#' Compute PCA
#' Function to compute PCA using
#' a modified algorithm combining the MM algorithm and Gradient Descent.
#' This algorithm is not using the full covariance matrix and instead decomposes
#' the data matrix to be memory efficient, even in the case of many variables.
#'
#' @param data_x An (nxp) data matrix
#' @param groups A (1xn) vector containing the groups of the observations
#' @param transformation The selected data transformation for the initial estimate.
#'  Default is "identity", this would correspond to using the sample covariance,
#'  has to be one of the following: "spearman", "wrapping"
#' @param loss_type The selected type of loss function. Default is "L2", corresponding to a least
#'  squares loss. The type has to be one of the following: "L2", "Huber", "Tukey", "LTS"
#' @param center logical: whether data should be centered before proceeding.
#'  Default is "TRUE", using the median
#'@param scale logical: whether data should be scaled before proceeding.
#'  Default is "TRUE", using the mad
#' @param alpha_x A positive real scalar containing a number between 0 and 1, indicating the elastic
#'  net parameter for data_x for each loading.
#' @param k highest order of maximum association you want to retrieve
#' @param tol desired accuracy for inner loop
#' @param epochs maximum number of iterations for inner loop
#' @param lr learning rate for gradient descent algorithm
#' @param lr_decay learning rate decay for gradient descent algorithm
#' @param criterion selection criterion for score function, default is Tradoff Product Optimization, TPO
#' @param penalties (optional) if given, no hyperparameter optimization is done,
#' but the given penalty parameters are used.
#' It has to be list structure, containing pen_x to be used for each loading
#' @param ... Additional parameters to be passed to the loss functions
#' @returns An object of class "PCA_result" that contains following entries:
#' @returns a A (pxk) vector of the estimated loadings
#' @returns measure A (1xk) vector of the explained variance
#' @returns phi A (nxk) matrix of the estimated scores corresponding a * data_x
#' @returns pen_x A (1xk) vector of the optimal penalty parameters for data_x
#' @returns alpha_x The value of the elastic net parameter alpha_x
#' @returns summary A summary of the hyperparameter optimization
#' @export
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats mad
#' @importFrom stats median
#' @importFrom stats rchisq
#' @importFrom stats sd
pcaSCRAMBLE <- function(data_x,
                        groups = NA,
                        transformation = "identity",
                        loss_type = "L2",
                        param = NA,
                        center = TRUE,
                        scale = TRUE,
                        alpha_x = NA,
                        k = NA,
                        tol = 1e-5,
                        lr = 1e-3,
                        epochs = 2000,
                        lr_decay = 1,
                        criterion = "TPO",
                        bounds = c(1e-3, 10),
                        penalties = NA) {

  p <- ncol(data_x)
  n <- nrow(data_x)

  if(is.na(k)){
    # per default, compute best possible rank k approximation
    # in case of really big matrices, it can speed up computation to provide a
    # smaller number for an initial result
    k <- min(p, n)
  }

  if(k > p){
    warning("k has to be <= min(p, n) ! k = min(p, n) is used")
    k <- min(p, n)
  }


  A <- matrix(NA, nrow = ncol(data_x), ncol = k)
  A_ortho <- matrix(NA, nrow = ncol(data_x), ncol = k)
  PHI <- matrix(NA, nrow = nrow(data_x), ncol = k)
  CORR <- rep(NA, k)
  CORR_ortho <- rep(NA, k)
  PEN_X <- NA

  SUMMARY <- list()

  wt <- rep(NA, n)
  # robscale <- NULL
  # robcenter <- NULL
  robcenter <- apply(data_x, 2, median, na.rm = TRUE)
  robscale <- apply(data_x, 2, robustbase::Qn, na.rm = TRUE)

  #data_og <- data_x

  # first, the data is standardized
  if (is.na(groups[1])){ # no groups
    if (isTRUE(center)){
      data_x <- scale(data_x,
                      center = robcenter,
                      scale = FALSE)
    }
    if (isTRUE(scale)){
      data_x <- scale(data_x,
                      center = FALSE,
                      scale = ifelse(robscale == 0, 1, robscale))
    }

  } else { # groups
      for (group in unique(groups)) {
        data_g <- data_x[which(groups == group),]
        robcenter <- rep(0, ncol(data_x))
        robscale <- rep(1, ncol(data_x))

      if (isTRUE(center)){
        robcenter <- apply(data_g, 2, median, na.rm = TRUE)
        data_g <- scale(data_g,
                        center = robcenter,
                        scale = FALSE)
        data_x[which(groups == group),] <- data_g
      }

      if (isTRUE(scale)){
        robscale <- apply(data_g, 2, stats::mad, na.rm = TRUE)
        data_g <- scale(data_g,
                        center = FALSE,
                        scale = ifelse(robscale == 0, 1, robscale))
        data_x[which(groups == group),] <- data_g
      }
     }
  }


  if (transformation == "identity") {
    rsvd <- svd(data_x)
    rsvd$d <- rsvd$d

  } else if (transformation == "spearman") {
    t_x <- apply(data_x, 2 , rank, ties.method = "average")
    t_x <-scale(t_x,
                center = apply(t_x, 2, mean),
                scale = apply(t_x, 2, sd))
    if (!isTRUE(scale)){
      t_x <- sweep(t_x, MARGIN=2, robscale, `*`)
    }

    if (!isTRUE(center)){
      t_x <- sweep(t_x, MARGIN=2, robcenter, `+`)
    }
    rsvd <- svd(t_x)
    rsvd$d <- rsvd$d

  }else if (transformation == "wrapping") {
    t_x <- cellWise::wrap(data_x,
                          locX = robcenter,
                          scaleX = robscale)$Xw
    rsvd <- svd(t_x)
    rsvd$d <- rsvd$d

  }
  # rsvd$v <- as.matrix(rsvd$v[,1:k])
  res <- apply(data_x %*% rsvd$v, 2, robustbase::Qn)
  order_res <- order(res, decreasing = TRUE)
  if (k > 1){
    rsvd$v <- rsvd$v[, order_res][, 1:k]
  } else {rsvd$v <- as.matrix(rsvd$v[, order_res][,1])}


    if (any(is.na(penalties))) {
      res_param <- bayesian_optimization_PCA_SVD(data_x, groups,
                                                 rsvd, n, p, k,
        alpha_x,
        loss_type, param,
        tol, lr, epochs,
        lr_decay = lr_decay,
        bounds_input = bounds,
        criterion = criterion
      )

      PEN_X <- as.numeric(res_param$best_params)

      SUMMARY <- res_param$summary
      # GAUPRO <- res_param$proc
    } else {
      PEN_X <- penalties$pen_x
    }
    res <- train_PCA_SVD(data_x,
                         groups,
                         rsvd,
                         n,
                         p,
                         PEN_X,
                         alpha_x,
                         loss_type,
                         param,
                         tol,
                         lr,
                         epochs,
                         lr_decay,
                         show_warnings = TRUE
                         )

    A <- res$best_v
    A_ortho <- res$best_v_ortho

    U <- res$best_u
    S <- res$best_s
    LOSS <- res$plot_loss

    #tot_var <- sum(apply(data_x, 2, robustbase::Qn))

    CORR <- res$best_measure
    CORR_ortho <- res$best_measure_ortho

    A <- A[,order(CORR, decreasing = TRUE)]
    A_ortho <- A_ortho[,order(CORR_ortho, decreasing = TRUE)]
    PHI <- as.matrix(data_x) %*% as.matrix(A)


    # RHO <- rep(NA, k)
    #
    # for (j in 1:k){
    #   M <- as.matrix(data_x) - as.matrix(data_x) %*% as.matrix(A[,order(CORR, decreasing = TRUE)][,1:j]) %*% t(as.matrix(A[,order(CORR, decreasing = TRUE)][, 1:j]))
    #   M_scale <- apply(M, 2, median)
    #   #M_center <- apply(M, 2, median)
    #   if (loss_type == "L2"){
    #     RHO[j] <- sum(M^2)
    #   } else if (loss_type == "LTS"){
    #     least_squares <- M^2
    #     threshold <- apply(least_squares, 2, quantile, probs = h)
    #     lts_loss <- ifelse(least_squares < threshold, least_squares, 0)
    #     loss <- sum(lts_loss)
    #     RHO[j] <- loss
    #   } else if (loss_type == "Huber"){
    #     M <- scale(M,
    #                center = FALSE,
    #                scale = M_scale)
    #     d <- 1
    #     huber <- d * sign(M) *M
    #     huber_loss <- ifelse(abs(M) < d, M^2, huber)
    #     loss <- sum(huber_loss)
    #     RHO[j] <-loss
    #   } else if (loss_type == "Tukey"){
    #     M <- scale(M,
    #                center = FALSE,
    #                scale = M_scale)
    #     k <- 1
    #     tukey <- (M/k)^2 * (3 - 3*(M/k)^2 + (M/k)^4)
    #     loss <- ifelse(abs(M) <= k, tukey, 1)
    #     RHO[j] <- sum(loss)
    #   }}


  return(structure(list(
    a = A,
    a_ortho = A_ortho,
    u = U,
    s = S,
    loss = LOSS,
    measure = as.numeric(CORR)[order(CORR, decreasing = TRUE)],
    measure_ortho = as.numeric(CORR_ortho)[order(CORR_ortho, decreasing = TRUE)],
    #explained_var = as.numeric(CORR)[order(CORR, decreasing = TRUE)] / tot_var,
    explained_var = cumsum(as.numeric(CORR_ortho)[order(CORR_ortho, decreasing = TRUE)]) / sum(as.numeric(CORR_ortho)[order(CORR_ortho, decreasing = TRUE)]),
    phi = PHI,
    pen_x = PEN_X[order(CORR, decreasing = TRUE)],
    summary = SUMMARY,
    #gauproc = GAUPRO,
    wt = wt,
    robscale = robscale,
    robcenter = robcenter,
    alpha_x = alpha_x[order(CORR, decreasing = TRUE)]
  ), class = "PCA_result"))
}

#' Plot PCA
#' Function to plot summary of PCA computation and hyperparameter optimization
#' @param x object of class PCA_result
#' @param ... additional parameters to be passed to ggplot2
#' @export
plot.PCA_result <- function(x, ...) {
  par(mfrow=c(1,2))
  # explained variance
  ev <- round(x$explained_var, 2)
  barplot(ev, xlab = "# PC", ylab = "Proportion of explained variance")
  plot(x$measure_ortho, xlab = "# PC", ylab = "Variances", type = "b", lty = 2, pch = 1)
  plot(x$phi[,1], x$phi[,2], xlab = "PC 1", ylab = "PC2", type = "p", pch = 1)
  diagplot(...,x)
  par(mfrow=c(1,3))
  plot(x$loss, type = "l", xlab = "Epoch", ylab = "Loss")
  if(length(x$summary) != 0){
    plot(sort(x$summary$pen_x), x$summary$Pred[order(x$summary$pen_x)], type = "l", xlab = "pen_x", ylab = "Total variance")
    abline(v = x$pen_x[1], lty = 2, col = "red")
    plot(sort(x$summary$pen_x), x$summary$Sparse_x[order(x$summary$pen_x)], type = "l", xlab = "pen_x", ylab = "Sparsity")
    abline(v = x$pen_x[1], lty = 2, col = "red")
  }
  par(mfrow=c(1,1))
  }

