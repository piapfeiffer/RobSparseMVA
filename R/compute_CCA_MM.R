#' Compute CCA
#' Function to compute CCA using combination of MM algorithm and Gradient Descent.
#'
#' @param data_x An (nxp) data matrix
#' @param data_y An (nxq) data matrix
#' @param method The selected method for computing the covariance.
#'  Default is "Pearson", referring to the sample covariance, has to be one
#'  of the following: "Spearman", "Kendall", "MCD", "MRCD", "OGK", "pairhuber",
#'  "quadrant", "Ledoit-Wolf".
#'  When a rank-based measure is selected, the covariance is computed as
#'  D . Cor . D, where D is a diagonal matrix containing the column-wise mad.
#' @param ... Additional parameters to be passed to the covariance estimators
#' @param nearPD logical, use Matrix::nearPD() on estimated covariance matrix
#' @param alpha_x (1xk) A positive real vector containing numbers between 0 (Ridge penalty) and 1 (LASSO penalty), indicating the elastic
#'  net parameter for data_x for association of order k.
#' @param alpha_y (1xk) A positive real vector containing numbers between 0 (Ridge penalty) and 1 (LASSO penalty), indicating the elastic
#'  net parameter for data_y for association of order k.
#' @param k highest order of maximum association you want to retrieve
#' @param tol desired accuracy for inner loop
#' @param epochs maximum number of iterations for inner loop
#' @param lr learning rate for gradient descent algorithm
#' @param lr_decay learning rate decay for gradient descent algorithm
#' @param criterion selection criterion for score function, default is Tradoff Product Optimization, TPO
#' @param penalties (optional) if given, no hyperparameter optimization is done,
#' but the given penalty parameters are used.
#' They have to be of list structure, containing pen_x and pen_y
#' @returns An object of class "CCA_result" that contains following entries:
#' @returns a A (pxk) vector of the estimated linear combinations
#'  corresponding to data_x.
#' @returns b A (qxk) vector of the estimated linear combinations
#'  corresponding to data_y.
#' @returns `measure` A (1xk) vector of the estimated maximum associations
#' @returns phi A (nxk) matrix of the estimated projections corresponding phi * data_x
#' @returns eta A (nxk) matrix of the estimated projections corresponding eta * data_y
#' @returns pen_x A (1xk) vector of the optimal penalty parameters for data_x
#' @returns pen_y A (1xk) vector of the optimal penalty parameters for data_y
#' @returns alpha_x The value of the elastic net parameter alpha_x
#' @returns alpha_y The value of the elastic net parameter alpha_y
#' @returns summary A summary of the hyperparameter optimization
#' @returns `method` The method used for computing the covariance matrix.
#' @export
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats mad
#' @importFrom stats median
#' @importFrom stats rchisq
#' @importFrom stats sd
ccaMM <- function(data_x, data_y,
                  method = "Pearson",
                  ...,
                  nearPD = FALSE,
                  alpha_x = NA,
                  alpha_y = NA,
                  k = 1,
                  tol = 1e-5,
                  lr = 1e-3,
                  epochs = 2000,
                  lr_decay = 1,
                  criterion = "TPO",
                  penalties = NA) {
  if (nrow(data_x) != nrow(data_y)) rlang::abort("Dimensions of x and y do not match", class = "data_error")
  if (anyNA(data_x) | anyNA(data_y)) rlang::abort("Data contains NA", class = "data_error")

  A <- matrix(NA, nrow = ncol(data_x), ncol = k)
  B <- matrix(NA, nrow = ncol(data_y), ncol = k)
  U <- matrix(NA, nrow = 2 * k + 2, ncol = k)

  PHI <- matrix(NA, nrow = nrow(data_x), ncol = k)
  ETA <- matrix(NA, nrow = nrow(data_y), ncol = k)

  CORR <- rep(NA, k)
  PEN_X <- rep(NA, k)
  PEN_Y <- rep(NA, k)
  SUMMARY <- list()

  p <- ncol(data_x)
  q <- ncol(data_y)
  n <- nrow(data_x)

  wt <- rep(NA, n)

  if (method == "Pearson") {
    C <- cov(cbind(data_x, data_y))
  } else if (method == "Spearman") {
    C <- sin(pi / 6 * cor(cbind(data_x, data_y), method = "spearman"))
    D <- diag(apply(cbind(data_x, data_y), 2, stats::mad))
    C <- D %*% C %*% D
    if (isTRUE(nearPD)) {
      C <- as.matrix(Matrix::nearPD(C)$mat)
    }
  } else if (method == "Kendall") {
    C <- sin(pi / 2 * cor(cbind(data_x, data_y), method = "kendall"))
    D <- diag(apply(cbind(data_x, data_y), 2, stats::mad))
    C <- D %*% C %*% D
    if (isTRUE(nearPD)) {
      C <- as.matrix(Matrix::nearPD(C)$mat)
    }
  } else if (method == "MCD") {
    cov_MCD <- rrcov::CovMcd(cbind(data_x, data_y), ...)
    C <- cov_MCD$cov
    wt <- cov_MCD$wt
  } else if (method == "MRCD") {
    cov_MRCD <- rrcov::CovMrcd(cbind(data_x, data_y), ...)
    C <- cov_MRCD$cov
    wt <- cov_MRCD$wt
  } else if (method == "OGK") {
    cov_OGK <- robustbase::covOGK(cbind(data_x, data_y), sigmamu = robustbase::s_mad)
    C <- cov_OGK$cov
    wt <- cov_OGK$weights
  } else if (method == "pairhuber") {
    C <- GSE::getScatter(GSE::HuberPairwise(cbind(data_x, data_y), psi=c("huber"), ...))
  } else if (method == "quadrant") {
    C <- GSE::getScatter(GSE::HuberPairwise(cbind(data_x, data_y), psi=c("sign"),...))
  } else if (method == "Ledoit-Wolf") {
    C <- BurStFin::var.shrink.eqcor(cbind(data_x, data_y), ...)
  }

  for (i in 1:k) {
    if (is.na(alpha_x[k])){
      alpha_x[k] <- 0
    }
    if (is.na(alpha_y[k])){
      alpha_y[k] <- 0
    }

    if (all(c(alpha_x == 0 & alpha_y == 0))){
      penalties <- list(pen_x = rep(sqrt(ncol(data_x)), k),
                        pen_y = rep(sqrt(ncol(data_y)), k))
    }
    if (any(is.na(penalties))) {
      res_param <- bayesian_optimization_CCA(C, p, q, n,
        alpha_x,
        alpha_y,
        i,
        low_a = as.matrix(A[, 1:max(1, (i - 1))]),
        low_b = as.matrix(B[, 1:max(1, (i - 1))]),
        tol, lr, epochs,
        lr_decay = lr_decay,
        criterion = criterion
      )

      PEN_X[i] <- res_param$best_params$pen_x
      PEN_Y[i] <- res_param$best_params$pen_y

      SUMMARY[[i]] <- res_param$summary
    } else {
      PEN_X[i] <- penalties$pen_x[i]
      PEN_Y[i] <- penalties$pen_y[i]
    }
    res <- train_CCA(C, p, q,
      PEN_X[i], PEN_Y[i], i,
      alpha_x[i], alpha_y[i],
      tol, lr, epochs,
      low_a = as.matrix(A[, 1:max(1, (i - 1))]),
      low_b = as.matrix(B[, 1:max(1, (i - 1))]),
      lr_decay,
      show_warnings = TRUE
    )

    A[, i] <- res$best_a
    B[, i] <- res$best_b
    U[1:(2 * i + 2), i] <- res$best_u

    PHI[, i] <- as.matrix(data_x) %*% as.matrix(res$best_a)
    ETA[, i] <- as.matrix(data_y) %*% as.matrix(res$best_b)

    CORR[i] <- res$best_measure
  }
  return(structure(list(
    a = A, b = B, measure = as.numeric(CORR),
    u = U,
    phi = PHI, eta = ETA,
    pen_x = PEN_X, pen_y = PEN_Y,
    summary = SUMMARY,
    method = method,
    cov = C,
    wt = wt,
    alpha_x = alpha_x, alpha_y = alpha_y
  ), class = "CCA_result"))
}

#' Plot CCA
#' Function to plot summary of CCA computation and hyperparameter optimization
#' @param x object of class CCA_result
#' @param ... additional parameters to be passed to ggplot2
#' @export
plot.CCA_result <- function(x, ...) {
  p_select <- list()
  p_sparse <- list()

  p_proj <- list()
  pen_x <- NULL
  pen_y <- NULL
  Pred <- NULL
  Score <- NULL
  Sparse_x <- NULL
  Sparse_y <- NULL
  phi <- NULL
  eta <- NULL

  for (i in 1:length(x$summary)) {
    p_select[[i]] <- ggplot2::ggplot(x$summary[[i]], ggplot2::aes(x = pen_x, y = pen_y, col = Pred, size = Score)) +
      ggplot2::geom_point() +
      ggplot2::geom_point(ggplot2::aes(x = x$pen_x[i], y = x$pen_y[i]), colour = "black", size = 10, shape = 1) +
      ggplot2::geom_point(ggplot2::aes(
        x = x$summary[[i]]$pen_x[which.max(x$summary[[i]]$Pred)],
        y = x$summary[[i]]$pen_y[which.max(x$summary[[i]]$Pred)]
      ), colour = "black", size = 10, shape = 0)

    xmin <- min(x$summary[[i]]$Sparse_x + x$summary[[i]]$Sparse_y)
    xmax <- x$summary[[i]]$Sparse_x[which.max(x$summary[[i]]$Score)] +
      x$summary[[i]]$Sparse_y[which.max(x$summary[[i]]$Score)]

    p_sparse[[i]] <- ggplot2::ggplot(x$summary[[i]], ggplot2::aes(x = Sparse_x + Sparse_y, y = Pred)) +
      ggplot2::geom_step() +
      ggplot2::geom_vline(ggplot2::aes(xintercept = xmax), linetype = 2, col = "blue") +
      ggplot2::ylim(0, 1)

    p_proj[[i]] <- ggplot2::ggplot(data.frame(phi = x$phi[, i], eta = x$eta[, i]), ggplot2::aes(x = phi, y = eta)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = x$measure[i]), col = "blue") +
      ggplot2::annotate("text",
        y = max(x$eta[, i]), x = min(x$phi[, i]),
        label = paste("  rho =", round(x$measure[i], 2)),
        col = "blue"
      )

    print(gridExtra::grid.arrange(p_select[[i]], p_sparse[[i]], p_proj[[i]],
      top = grid::textGrob(paste("order:", i), gp = grid::gpar(fontsize = 20)), ncol = 3
    ))
  }
}
