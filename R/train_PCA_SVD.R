#' MM algorithm for PCA objects
#' implements MM algorithms for robust and sparse PCA
#' @param C estimated covariance matrix
#' @param p number of variables in data_x
#' @param pen_x (1xk) A positive real vector containing the penalty parameter for
#'  data_x for association of order k.
#' @param alpha_x (1xk) A positive real vector containing numbers between 0 and
#' 1, indicating the elastic
#'  net parameter for data_x for association of order k.
#' @param tol tolerance for training algorithm
#' @param lr  learning rate for gradient descent training algorithm
#' @param epochs maximum number of epochs for gradient descent training algorithm
#' @param low_a A (px(k-1)) matrix of lower order linear combinations
#' @param lr_decay parameter for
#' @param show_warnings logical
#' @param init_ortho logical
#' @export
#' @importFrom utils head
train_PCA_SVD <- function(X,
                          groups,
                          rsvd, n, p,
                          pen_x,
                          alpha_x,
                          loss_type,
                          h,
                          tol, lr, epochs,
                          lr_decay,
                          show_warnings = FALSE,
                          init_ortho = TRUE) {

  best_u <- NA
  best_v <- NA
  best_s <- NA
  best_m <- NA
  old_best_u <- NA
  old_best_v <- NA
  old_best_s <- NA

  plot_loss <- rep(NA, epochs)

    model <- model_PCA_SVD(X,
                           groups,
                           rsvd, n, p,
                       pen_x = pen_x,
                       alpha_x = alpha_x,
                       loss_type = loss_type,
                       h = h,
                       warm_v = best_v,
                       warm_m = best_m)

    best_v <- as.matrix(model$parameters$V$detach())

    best_loss <- 1e16
    best_measure <- 0
    change <- 1e16
    count_ep <- 0

    old_v <- best_v
    old_loss <- best_loss

    lr <- lr
    optimizer2 <- optim_sgd_manifold(model$parameters$V, lr = lr, manifold = TRUE)
    lmbda <- function(epoch) lr_decay
    scheduler2 <- torch::lr_multiplicative(optimizer2, lr_lambda = lmbda)

    model$zero_grad()

    diff_log_v <- rstackdeque::rdeque()

    while (as.numeric(change) > tol && count_ep < epochs) {
      optimizer2$zero_grad()

      loss <- model(X)
      plot_loss[count_ep] <- as.numeric(loss)

      if (as.numeric(loss) < best_loss) {
        best_loss <- as.numeric(loss)
        best_v <- as.matrix(model$parameters$V$detach())

        best_measure <- model$evaluate(X, best_v)
      }

      loss$backward()

      optimizer2$step()
      scheduler2$step()

      check <- t(best_v) %*% best_v



      if(norm(old_v, "1") != 0){
        diff_log_v <- rstackdeque::insert_front(diff_log_v,
                                                norm(as.matrix(model$parameters$V$detach()) - old_v, "F") /
                                                  norm(old_v, "F"))
      } else {
        diff_log_v <- rstackdeque::insert_front(diff_log_v, 0)
      }

      old_v <- as.matrix(model$parameters$V$detach())

      count_ep <- count_ep + 1
      change <- abs(as.numeric(old_loss - loss))
      old_loss <- loss
      #print(loss)


    } # end inner loop


    norm_vec <- function(x) sqrt(sum(x^2))

    threshold_v <- mean(as.numeric(as.list(head(diff_log_v, 10)))) + 2 * sd(as.numeric(as.list(head(diff_log_v, 10))))
    #print(threshold_v)
    best_v_ortho <- best_v
    best_measure_ortho <- best_measure
    best_v_thresh <- ifelse(abs(best_v) > threshold_v, best_v, 0)
       if (norm(as.matrix(best_v_thresh), "1") > 0){
        best_v <- best_v_thresh

        scale <- ifelse(apply(best_v, 2, norm_vec)==0, 1, apply(best_v, 2, norm_vec))

        best_v <- scale(best_v,
                        center = FALSE,
                        scale = scale)
        best_measure <- model$evaluate(X, best_v)
       }


    if(as.numeric(change) > tol){
        if (isTRUE(show_warnings)){
        warning(paste("Loop did not converge. Try more iterations. Change:", as.numeric(change)))
        } else {
          best_measure <- 0
        }
  }

  return(list(best_measure = best_measure,
              best_loss = best_loss,
              best_v = best_v,
              best_v_ortho = best_v_ortho,
              plot_loss = plot_loss,
              best_measure_ortho = best_measure_ortho
  ))
}
