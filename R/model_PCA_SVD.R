#' NN module for robust and sparse PCA via SVD
#' @param X data input
#' @param rsvd object containing SVD decompositin of R
#' @param n number of
#' @param p number of variables in data_x
#' @param pen_x (1xk) A positive real vector containing the penalty parameter for
#'  data_x for association of order k.
#' @param order The order of association to be computed
#' @param alpha_x (1xk) A positive real vector containing numbers between 0 and
#' 1, indicating the elastic
#'  net parameter for data_x for association of order k.
#' @param warm_u A vector of the currently best linear combination
#' @param warm_v A vector of the currently best linear combination
#' @param warm_s A  scalar of the currently best singular value
#' @param warm_m A vector of the currently best multiplier values
#' @param low_u A (px(k-1)) matrix of lower order linear combinations

model_PCA_SVD <- torch::nn_module(
  initialize = function(X,
                        groups,
                        rsvd,
                        n,
                        p,
                        pen_x,
                        alpha_x,
                        loss_type,
                        h,
                        warm_v = NA,
                        warm_m = NA) {
    self$loss_type <- loss_type
    self$p <- p
    self$n <- n
    self$h <- h
    self$pen_x <- pen_x #TODO: check whether dimensions are correct
    self$groups <- groups #TODO: check whether dimensions are correct

    self$alpha_x <- alpha_x

    if (is.na(warm_v[1])) {
        v <- rsvd$v
        self$V <- torch::nn_parameter(torch::torch_tensor(v,
                                                          dtype = torch::torch_double()), requires_grad = TRUE)
        # self$V <- torch::nn_parameter(torch::torch_eye(nrow(v), ncol(v),
        #                                                dtype = torch::torch_double()),
        #                               requires_grad = TRUE)

      } else {
        v <- warm_v
        self$V <- torch::nn_parameter(torch::torch_tensor(v,
                                                          dtype = torch::torch_double()), requires_grad = TRUE)
      }
    self$TV <- apply(X, 2, robustbase::Qn)
    # for (j in 1:ncol(v)){
    #   self$TV[j] <- sum(apply(X %*% v[,1:(j-1)], 2, robustbase::Qn))
    # }

  }, forward = function(X) {

    m1 <- torch::torch_tensor(rep(0.0, ncol(self$V)),
                              dtype = torch::torch_double())
    for (j in 1:ncol(self$V)){
      m1[j] <- penalty_fun_unc(self$V[,j], self$pen_x * self$TV[j], self$alpha_x)
    }

    G <- torch::torch_cat(c(torch::torch_unsqueeze(m1, 2)))

    B <-  torch::torch_matmul(torch::torch_tensor(X, dtype = torch::torch_double()),
                              torch::torch_matmul(self$V,
                                                  torch::torch_transpose(self$V, 1, -1)))
    if (self$loss_type == "L2") {
      loss <- torch::torch_mean(torch::torch_square(X - B)) + torch::torch_sum(G)

    } else if (self$loss_type == "Huber") {
      loss <- loss_huber(X - B) + torch::torch_sum(G)

    } else if (self$loss_type == "Tukey") {
      loss <- loss_tukey(X - B) + torch::torch_sum(G)

    } else if (self$loss_type == "LTS") {
      if (is.na(self$groups[1])){
        loss <- loss_lts(X - B, self$h) + torch::torch_sum(G)
      }
      else{
        loss <- torch::torch_sum(G)
        for (group in unique(self$groups)){
          XG <- X[which(self$groups == group),]
          BG <-  torch::torch_matmul(torch::torch_tensor(XG, dtype = torch::torch_double()),
                                    torch::torch_matmul(self$V,
                                                        torch::torch_transpose(self$V, 1, -1)))
          loss <- loss + loss_lts(XG - BG, self$h)
        }
      }

    }else {
      warning(paste("Loss type:",self$loss_type, "not implemented, standard L2 loss is used instead!"))
      loss <- torch::torch_square(torch::torch_norm(X - B)) + torch::torch_sum(G)
    }

     loss
  }, evaluate = function(X, V) {
    #res <- apply(X %*% V, 2, mad)
    res <- apply(X %*% V, 2, robustbase::Qn)
    return(res)
  }
)
