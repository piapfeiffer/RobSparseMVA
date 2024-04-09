#' Get principal angles.
#'
#' Calculate minimal angle between subspace a and b
#' @param a  A matrix.
#' @param b  A matrix.
#' @returns  A list of angles.
#' @export
principal_angles <- function(a, b) {
  angles <- matrix(0, ncol = ncol(a), nrow = 1)
  qa <- qr.Q(qr(a))
  qb <- qr.Q(qr(b))
  C <- svd(t(qa) %*% qb)$d
  rkA <- qr(a)$rank
  rkB <- qr(b)$rank
  if (rkA <= rkB) {
    B <- qb - qa %*% (t(qa) %*% qb)
  } else {
    B <- qa - qb %*% (t(qb) %*% qa)
  }
  S <- svd(B)$d
  S <- sort(S)

  if (min(rkA, rkB) > 0) {
    for (i in 1:min(rkA, rkB)) {
      if (C[i]^2 < 0.5) {
        angles[1, i] <- acos(C[i])
      } else if (S[i]^2 <= 0.5) {
        angles[1, i] <- asin(S[i])
      }
    }
    angles <- t(angles)
  }
  out <- list(angles = angles)
}


#' True positive rate
#'
#' This function compares the structure of two matrices A and A_true.
#'
#' @param A_true A matrix of the ground truth.
#' @param A A matrix.
#' @returns The number of correctly identified non-zero components.
#' @export
TPR <- function(A_true, A) {
  out <- sum(A_true != 0 & A != 0) / sum(A_true != 0)
}


#' True negative rate
#'
#' This function compares the structure of two matrices A and A_true.
#'
#' @param A_true   A matrix of the ground truth.
#' @param A        A matrix.
#' @returns The number of correctly identified zero components.
#' @export
TNR <- function(A_true, A) {
  out <- sum(A_true == 0 & A == 0) / sum(A_true == 0)
}

#' Higher order constraint for robust and sparse PCA / CCA
#'
#' Function to compute higher order constraint for a tensor
#'
#' @param z       A (1xp) torch::torch_tensor
#' @param low_z   A (kxp) torch::torch_tensor
#' @param C       A (pxp) matrix
higher_order <- function(z, low_z, C = NULL) {

  if (is.null(C)) {
    low_u <- -torch::torch_sum(torch::torch_abs(torch::torch_matmul(torch::torch_transpose(z, 1, -1), low_z)), 2)
  } else{
    # low_u <- torch::torch_matmul(torch::torch_transpose(z, 1, -1), torch::torch_matmul(C, low_z))
    # return(torch::torch_unsqueeze(torch::torch_flatten(torch::torch_transpose(low_u, 1, 2)),2))
    low_u <- -torch::torch_sum(torch::torch_abs(torch::torch_matmul(torch::torch_transpose(z, 1, -1), torch::torch_matmul(C, low_z))), 2)}
  # print(torch::torch_unsqueeze(low_u, 2))
  return(torch::torch_unsqueeze(low_u, 2))
}

#' Domain constraint for robust and sparse PCA / CCA
#'
#' Function to compute domain constraint for a tensor
#'
#' @param z    A (1xp) torch::torch_tensor
#' @param C    A (pxp) quadratic matrix
#' @returns z_constrained, the value of the equality constraints
constraint_fun <- function(z, C = NULL) {
  if (is.null(C)) {
    z_constrained <- torch::torch_matmul(torch::torch_transpose(z, 1, -1), z) - 1
  } else {
    z_constrained <- torch::torch_matmul(torch::torch_transpose(z, 1, -1), torch::torch_matmul(C, z)) - 1
  }
  if (as.numeric(z_constrained) > 0) {
    return(z_constrained)
  } else {
    return(torch::torch_unsqueeze(torch::torch_tensor(c(0)), 2))
  }
}

#' Penalty constraints for robust and sparse PCA / CCA
#'
#' Function to compute elastic net penalty for a tensor
#'
#' @param z             A torch::torch_tensor
#' @param pen           penalty parameter
#' @param alpha         elastic net mixing parameter
#' @returns              A 1x1 torch::torch_tensor
penalty_fun <- function(z, pen, alpha) {
  l1 <- torch::torch_sum(torch::torch_abs(z))
  l2 <- torch::torch_square(torch::torch_norm(z))
  penalty <- alpha * l1 + (1 - alpha) * l2 - pen

  if (as.numeric(penalty) > 0) {
    return(torch::torch_unsqueeze(penalty, 2))
  } else {
    return(torch::torch_unsqueeze(torch::torch_tensor(c(0)), 2))
  }
}

#' Penalty constraints for robust and sparse PCA via SVD
#'
#' Function to compute elastic net penalty for a tensor
#'
#' @param z             A torch::torch_tensor
#' @param pen           penalty parameter
#' @param alpha         elastic net mixing parameter
#' @returns              A 1x1 torch::torch_tensor
penalty_fun_unc <- function(z, pen, alpha) {
  l1 <- torch::torch_sum(torch::torch_tanh(1000*z) * z)
  l2 <- torch::torch_square(torch::torch_norm(z))
  penalty <- (alpha * l1 + (1 - alpha) * l2) * pen

  return(torch::torch_unsqueeze(penalty, 2))

}

#' Loss function for robust and sparse PCA via SVD
#'
#' Function to compute Pseudo-Huber loss for a given matrix.
#'
#' @param M             A matrix torch::torch_tensor
#' @param d             Scalar, parameter for Huber loss
#' @returns             A 1x1 torch::torch_tensor
loss_huber <- function(M, d = 1.35) {
  # M <- scale(M,
  #            center = FALSE,
  #            scale = apply(abs(M), 2, median))
  sigma <- torch::torch_median(torch::torch_abs(M), 1, keepdim = TRUE)[[1]]
  M <- M / sigma
  #least_squares <- M^2
  #huber <- torch::torch_tanh(1000*M) * M * d
  #loss <- torch::torch_sum(torch::torch_where(torch::torch_abs(M) <= d, least_squares, huber))
  loss <- torch::torch_mean(d^2 * (torch::torch_sqrt(1 + (M/d)^2) - 1) * sigma^2)
  return(loss)
}

#' Loss function for robust and sparse PCA via SVD
#'
#' Function to compute Tukey loss for a given matrix.
#'
#' @param M             A matrix torch::torch_tensor
#' @param k             Scalar, parameter for Huber loss
#' @returns             A 1x1 torch::torch_tensor
loss_tukey <- function(M, k = 1.35) {
  # M <- scale(M,
  #            center = FALSE,
  #            scale = apply(abs(M), 2, median))
  sigma <- torch::torch_median(torch::torch_abs(M), 1, keepdim = TRUE)[[1]]
  M <- M / sigma
  tukey <- (M/k)^2 * (3 - 3*(M/k)^2 + (M/k)^4)
  loss <- torch::torch_mean(torch::torch_where(torch::torch_abs(M) <= k, tukey, 1.) * sigma^2)
  return(loss)
}

#' Loss function for robust and sparse PCA via SVD
#'
#' Function to compute trimmes loss for a given matrix.
#'
#' @param M             A matrix torch::torch_tensor
#' @param h             Trimming parameter (robustness) of observations to use per column
#' @returns             A 1x1 torch::torch_tensor
loss_lts <- function(M, h = 0.75) {
  least_squares <- M^2
  threshold <- torch::torch_quantile(M,
                                     torch::torch_tensor(h, dtype = torch::torch_double()),
                                     dim = 1,
                                     keepdim = TRUE,
                                     interpolation = "nearest")
  lts_loss <- torch::torch_where(least_squares < threshold, least_squares, 0)
  loss <- torch::torch_mean(lts_loss)
  return(loss)
}


#' F1 score
#'
#' Function to compute F1 score for outlier detection, given the true positions
#'
#' @param ind_outl     position of outliers
#' @param flag_outl    flagged cells
f1_score <- function(ind_outl, flag_outl) {
  TP <- length(intersect(flag_outl, ind_outl)) #TP
  FP <- length(setdiff(flag_outl, ind_outl)) #FP
  FN <- length(setdiff(ind_outl, flag_outl)) #FN

  f1_score <- TP / (TP + 0.5 * (FP + FN))
}
