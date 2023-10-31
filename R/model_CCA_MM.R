#' NN module for robust and sparse CCA
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
#' @param rho A real number. The current value of the penalty strength in the
#' Augmented Lagrangian function.
#' @param warm_a A (px1) vector of the currently best linear combination
#' @param warm_b A (qx1) vector of the currently best linear combination
#' @param warm_u A vector of the currently best multiplier values
#' @param low_a A (px(k-1)) matrix of lower order linear combinations
#' @param low_b A (qx(k-1)) matrix of lower order linear combinations
#' @param init_ortho logical: whether initialization should be orthogonal
#' for higher-order directions
model_CCA <- torch::nn_module(
  initialize = function(C, p, q,
                        pen_x, pen_y,
                        order,
                        alpha_x, alpha_y,
                        rho,
                        warm_a = NA, warm_b = NA, warm_u = NA,
                        low_a = NA, low_b = NA,
                        init_ortho = TRUE) {
    self$rho <- rho
    self$p <- p
    self$q <- q
    self$pen_x <- pen_x
    self$pen_y <- pen_y

    self$order <- order
    self$alpha_x <- alpha_x
    self$alpha_y <- alpha_y

    if (order > 1) {
      self$low_a <- torch::torch_tensor(low_a)
      self$low_b <- torch::torch_tensor(low_b)
    }

    if (p > 1) {
      if (is.na(warm_a[1])) {
        if (order > 1) {
          if (isTRUE(init_ortho)) {
            start_a <- apply(t(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)])), 2, mean)
            nullspace_a <- pracma::nullspace(t(C[1:self$p, 1:self$p] %*% low_a))
            nom_a <- 0
            denom_a <- 0
            proj_a <- rep(0, length(start_a))
            for (i in 1:ncol(nullspace_a)) {
              vec <- as.numeric(nullspace_a[, i])
              nom_a <- nom_a + as.numeric(start_a %*% vec)
              denom_a <- denom_a + as.numeric(t(vec) %*% vec)
              proj_a <- proj_a + (nom_a / denom_a) * vec
            }
            a <- proj_a
          } else {
            a <- apply(t(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)])), 2, mean)
          }

          denom_a <- sqrt(t(a) %*% as.matrix(C[1:self$p, 1:self$p]) %*% a)[1, 1]
          if (denom_a > 0) {
            a <- a / denom_a
          }
        } else {
          a <- apply(t(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)])), 2, mean)
          denom_a <- sqrt(t(a) %*% as.matrix(C[1:self$p, 1:self$p]) %*% a)[1, 1]
          if (denom_a > 0) {
            a <- a / denom_a
          }
        }
        self$a <- torch::nn_parameter(torch::torch_unsqueeze(torch::torch_tensor(a), 2), requires_grad = TRUE)
      } else {
        a <- warm_a
        self$a <- torch::nn_parameter(torch::torch_tensor(a), requires_grad = TRUE)
      }
    } else {
      self$a <- torch::nn_parameter(torch::torch_unsqueeze(torch::torch_tensor(1), 2), requires_grad = TRUE)
    }

    if (q > 1) {
      if (is.na(warm_b[1])) {
        if (order > 1) {
          if (isTRUE(init_ortho)) {
            start_b <- apply(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)]), 2, mean)
            nullspace_b <- pracma::nullspace(t(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)] %*% low_b))
            nom_b <- 0
            denom_b <- 0
            proj_b <- rep(0, length(start_b))
            for (i in 1:ncol(nullspace_b)) {
              vec <- as.numeric(nullspace_b[, i])
              nom_b <- nom_b + as.numeric(start_b %*% vec)
              denom_b <- denom_b + as.numeric(t(vec) %*% vec)
              proj_b <- proj_b + (nom_b / denom_b) * vec
            }
            b <- proj_b
          } else {
            b <- apply(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)]), 2, mean)
          }

          denom_b <- sqrt(t(b) %*% as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)]) %*% b)[1, 1]
          if (denom_b > 0) {
            b <- b / denom_b
          }
        } else {
          b <- apply(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)]), 2, mean)
          denom_b <- sqrt(t(b) %*% as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)]) %*% b)[1, 1]
          if (denom_b > 0) {
            b <- b / denom_b
          }
        }
        self$b <- torch::nn_parameter(torch::torch_unsqueeze(torch::torch_tensor(b), 2), requires_grad = TRUE)
      } else {
        b <- warm_b
        self$b <- torch::nn_parameter(torch::torch_tensor(b), requires_grad = TRUE)
      }
    } else {
      self$b <- torch::nn_parameter(torch::torch_unsqueeze(torch::torch_tensor(1), 2), requires_grad = TRUE)
    }

    current_corr <- as.numeric(torch::torch_matmul(
      torch::torch_transpose(self$a, 1, -1),
      torch::torch_matmul(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)]), self$b)
    ))
    if (current_corr < 0) {
      self$a <- -self$a
    }

    if (order > 1) {
      if (is.na(warm_u[1])) {
        u1_a <- constraint_fun(self$a$detach(), as.matrix(C[1:self$p, 1:self$p]))
        u1_b <- constraint_fun(self$b$detach(), as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)]))
        u2_a <- penalty_fun(self$a$detach(), self$pen_x, self$alpha_x)
        u2_b <- penalty_fun(self$b$detach(), self$pen_y, self$alpha_y)
        u3_a <- higher_order(self$a$detach(), self$low_a, torch::torch_tensor(as.matrix(C[1:self$p, 1:self$p])))
        u3_b <- higher_order(self$b$detach(), self$low_b, torch::torch_tensor(as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)])))
        self$u <- torch::torch_stack(c(u1_a, u1_b, u2_a, u2_b, u3_a, u3_b))
      } else {
        self$u <- torch::torch_tensor(warm_u)
      }
    } else {
      if (is.na(warm_u[1])) {
        u1_a <- constraint_fun(self$a$detach(), as.matrix(C[1:self$p, 1:self$p]))
        u1_b <- constraint_fun(self$b$detach(), as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)]))
        u2_a <- penalty_fun(self$a$detach(), self$pen_x, self$alpha_x)
        u2_b <- penalty_fun(self$b$detach(), self$pen_y, self$alpha_y)
        self$u <- torch::torch_stack(c(u1_a, u1_b, u2_a, u2_b))
      } else {
        self$u <- torch::torch_tensor(warm_u)
      }
    }
  }, forward = function(rho, C) {
    if (self$order > 1) {
      G <- torch::torch_stack(c(
        constraint_fun(
          self$a,
          as.matrix(C[1:self$p, 1:self$p])
        ),
        constraint_fun(
          self$b,
          as.matrix(C[
            (self$p + 1):(self$p + self$q),
            (self$p + 1):(self$p + self$q)
          ])
        ),
        penalty_fun(self$a, self$pen_x, self$alpha_x),
        penalty_fun(self$b, self$pen_y, self$alpha_y),
        higher_order(self$a, self$low_a, as.matrix(C[1:self$p, 1:self$p])),
        higher_order(self$b, self$low_b, as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)]))
      ))[, , 1]
    } else {
      G <- torch::torch_stack(c(
        constraint_fun(self$a, as.matrix(C[1:self$p, 1:self$p])),
        constraint_fun(self$b, as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)])),
        penalty_fun(self$a, self$pen_x, self$alpha_x),
        penalty_fun(self$b, self$pen_y, self$alpha_y)
      ))[, , 1]
    }
    -torch::torch_abs(torch::torch_matmul(
      torch::torch_transpose(self$a, 1, -1),
      torch::torch_matmul(as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)]), self$b)
    )) +
      torch::torch_matmul(torch::torch_transpose(self$u, 1, -1), G) +
      rho / 2 * torch::torch_square(torch::torch_norm(G))
  }, evaluate = function(a, b, C) {
    res <- (t(a) %*% as.matrix(C[1:self$p, (self$p + 1):(self$p + self$q)]) %*% b) /
      (sqrt(t(a) %*% as.matrix(C[1:self$p, 1:self$p]) %*% a) * sqrt(t(b) %*% as.matrix(C[(self$p + 1):(self$p + self$q), (self$p + 1):(self$p + self$q)]) %*% b))
    if (is.na(res) == TRUE) {
      return(0)
    } else {
      return(res)
    }
  }
)
