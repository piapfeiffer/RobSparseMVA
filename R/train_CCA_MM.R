#' MM algorithm for CCA objects
#' implements MM algorithms for robust and sparse CCA
#' @param C estimated covariance matrix
#' @param p number of variables in data_x
#' @param q number of variables in data_y
#' @param pen_x (1xk) A positive real vector containing the penalty parameter for
#'  data_x for association of order k.
#' @param pen_y (1xk) A positive real vector containing the penalty parameter for
#' data_y for association of order k.
#' @param order The order of association to be computed
#' @param alpha_x (1xk) A positive real vector containing numbers between 0 and
#' 1, indicating the elastic
#'  net parameter for data_x for association of order k.
#' @param alpha_y (1xk) A positive real vector containing numbers between 0 and
#' 1, indicating the elastic
#'  net parameter for data_y for association of order k.
#' @param tol tolerance for training algorithm
#' @param lr  learning rate for gradient descent training algorithm
#' @param epochs maximum number of epochs for gradient descent training algorithm
#' @param low_a A (px(k-1)) matrix of lower order linear combinations
#' @param low_b A (qx(k-1)) matrix of lower order linear combinations
#' @param lr_decay parameter for
#' @param show_warnings logical
#' @param init_ortho logical
#' @export
#' @importFrom utils head
train_CCA <- function(C, p, q,
                      pen_x, pen_y, order,
                      alpha_x, alpha_y,
                      tol, lr, epochs,
                      low_a, low_b,
                      lr_decay,
                      show_warnings = FALSE,
                      init_ortho = TRUE) {

  outer_count <- 0
  outer_change <- 1e16
  rho <- 1
  best_u <- NA
  best_a <- NA
  best_b <- NA
  old_best_a <- NA
  old_best_b <- NA
  outer_epochs <- 30
  loss_outer <- rep(NA, outer_epochs + 1)
  loss_outer[1] <- 1e16
  diff <- 1e16
  threshold_a <- 0
  threshold_b <- 0

  while (as.numeric(diff) > tol && outer_count < outer_epochs) {
    old_loss_outer <- loss_outer
    model <- model_CCA(C, p, q,
                       pen_x = pen_x,
                       pen_y = pen_y,
                       order = order,
                       alpha_x = alpha_x,
                       alpha_y = alpha_y,
                       rho = rho,
                       warm_a = best_a,
                       warm_b = best_b,
                       warm_u = best_u,
                       low_a = low_a,
                       low_b = low_b,
                       init_ortho = init_ortho)

    best_a <- as.matrix(model$parameters$a$detach())
    best_b <- as.matrix(model$parameters$b$detach())

    old_u <- model$u

    old_h <- old_u

    best_loss <- 1e16
    best_measure <- 0
    change <- 1e16
    count_ep <- 0

    old_a <- as.matrix(model$parameters$a)

    old_b <- as.matrix(model$parameters$b)

    old_loss <- best_loss

    lr <- lr
    optimizer <- torch::optim_adam(model$parameters, lr = lr, amsgrad = TRUE)
    lmbda <- function(epoch) lr_decay
    scheduler <- torch::lr_multiplicative(optimizer, lr_lambda = lmbda)

    model$zero_grad()

    diff_log_a <- rstackdeque::rdeque()
    diff_log_b <- rstackdeque::rdeque()

    while (as.numeric(change) > tol && count_ep < epochs) {
      optimizer$zero_grad()

      loss <- model(rho, C)
      if(is.na(as.numeric(loss))){
        print(as.numeric(loss))
        print(as.numeric(model$a))
        print(as.numeric(model$b))
      }


      if (as.numeric(loss) < best_loss) {
        best_loss <- as.numeric(loss)

        best_a <- as.matrix(model$parameters$a$detach())
        best_b <- as.matrix(model$parameters$b$detach())
        best_u <- as.matrix(model$u$detach())

        best_measure <- model$evaluate(best_a,
                                       best_b,
                                       C)
      }

      loss$backward()

      optimizer$step()
      scheduler$step()

      count_ep <- count_ep + 1

      change <- abs(as.numeric(old_loss - loss))

      if(norm(old_a, "1") != 0){
        diff_log_a <- rstackdeque::insert_front(diff_log_a,
                                                norm(as.matrix(model$parameters$a$detach()) - old_a, "F") /
                                                  norm(old_a, "F"))
      } else {
        diff_log_a <- rstackdeque::insert_front(diff_log_a, 0)
      }
      if (norm(old_b, "1") != 0){
        diff_log_b <- rstackdeque::insert_front(diff_log_b,
                                                norm(as.matrix(model$parameters$b$detach()) - old_b, "F") /
                                                  norm(old_b, "F"))
      } else {
        diff_log_b <- rstackdeque::insert_front(diff_log_b, 0)
      }

      old_a <- as.matrix(model$parameters$a$detach())
      old_b <- as.matrix(model$parameters$b$detach())

      old_loss <- loss

    } # end inner loop

      threshold_a <- mean(as.numeric(as.list(head(diff_log_a, 10)))) + 2 * sd(as.numeric(as.list(head(diff_log_a, 10))))
      if (norm(best_a, "F") > 0) {
        best_a_thresh <- ifelse(abs(best_a) > threshold_a, best_a, 0)
      }
      if (norm(best_a_thresh, "1") > 0){
        best_a <- best_a_thresh
      }

      threshold_b <- mean(as.numeric(as.list(head(diff_log_b, 10)))) + 2 * sd(as.numeric(as.list(head(diff_log_b, 10))))
      if (norm(best_b, "F") > 0) {
        best_b_thresh <- ifelse(abs(best_b) > threshold_b, best_b, 0)
      }
      if (norm(best_b_thresh, "1") > 0){
        best_b <- best_b_thresh
      }

    denom_best_a <- t(best_a) %*% C[1:p, 1:p] %*% best_a
    if (denom_best_a > 0) {
      best_a <- best_a / sqrt(denom_best_a[1,1])
    }

    denom_best_b <- t(best_b) %*% C[(p+1):(p + q), (p+1):(p+q)] %*% best_b
    if (denom_best_b > 0) {
      best_b <- best_b / sqrt(denom_best_b[1,1])
    }

    best_measure <- model$evaluate(best_a,
                                   best_b,
                                   C)

    torch::with_no_grad({
    if (order > 1){
      model$u <- model$u + rho * torch::torch_stack(c(constraint_fun(torch::torch_tensor(best_a), C[1:p, 1:p]),
                                                      constraint_fun(torch::torch_tensor(best_b), C[(p+1):(p + q), (p+1):(p+q)]),
                                                      penalty_fun(torch::torch_tensor(best_a),pen_x, alpha_x),
                                                      penalty_fun(torch::torch_tensor(best_b),pen_y, alpha_y),
                                                      higher_order(torch::torch_tensor(best_a),model$low_a, C[1:p, 1:p]),
                                                      higher_order(torch::torch_tensor(best_b),model$low_b, C[(p+1):(p + q), (p+1):(p+q)])))
    } else {
      model$u <- model$u + rho * torch::torch_stack(c(constraint_fun(torch::torch_tensor(best_a), C[1:p, 1:p]),
                                                      constraint_fun(torch::torch_tensor(best_b), C[(p+1):(p + q), (p+1):(p+q)]),
                                                      penalty_fun(torch::torch_tensor(best_a),pen_x, alpha_x),
                                                      penalty_fun(torch::torch_tensor(best_b),pen_y, alpha_y)))
    }
      new_u <- model$u
      new_h <- (new_u - old_u) / rho
    })

    outer_change <- as.numeric(torch::torch_norm((old_u - new_u)/rho))

    if (0.25 * as.numeric(torch::torch_norm(old_h, 1)) < as.numeric(torch::torch_norm(new_h, 1))){
      if (rho < 1e14){
        rho <- 10 * rho
      }
    }
    loss_outer[outer_count + 2] <- as.numeric(outer_change)
    diff <- loss_outer[outer_count + 1] - loss_outer[outer_count + 2]
    outer_count <- outer_count + 1
  }


    if(as.numeric(diff) > tol){
        if (isTRUE(show_warnings)){
        warning(paste("Outer Loop did not converge. Try more iterations. Change:", as.numeric(outer_change)))
        } else {
          best_measure <- 0
        }
    }
    if(as.numeric(change) > tol){
        if (isTRUE(show_warnings)){
        warning(paste("Inner Loop did not converge. Try more iterations. Change:", as.numeric(change)))
        } else {
          best_measure <- 0
        }
  }

  return(list(best_measure = best_measure,
              best_loss = best_loss,
              best_a = best_a,
              best_b = best_b,
              best_u = best_u
  ))
}
