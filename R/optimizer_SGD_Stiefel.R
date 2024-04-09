optim_sgd_manifold <- torch::optimizer(
  initialize = function(params, lr, manifold, threshold = NA) {
    defaults <- list(
      lr = lr,
      manifold = manifold,
      threshold = threshold
    )
    super$initialize(params, defaults)
  },
  step = function() {

    torch::with_no_grad({

        for (g in seq_along(self$param_groups)) {
          group <- self$param_groups[[g]]
          if (isTRUE(group$manifold)){

            for (p in seq_along(group$params)) {
              param <- group$params[[p]]

              if (is.null(param$grad) || torch::is_undefined_tensor(param$grad)) {
                next
              }

              manifold_grad <- -group$lr * (param$grad - torch::torch_matmul(param,
                                                                torch::torch_matmul(
                                                                  torch::torch_transpose(param$grad, 1, -1),
                                                                  param)
                                                                ))

              a1 <- torch::torch_matmul(manifold_grad,
                                        torch::torch_transpose(param, 1, -1)
                                        )
              a2 <-  torch::torch_matmul(param,
                                         torch::torch_transpose(manifold_grad, 1, -1))
              a <- a1 - a2
              rhs <- param + 1/2 * torch::torch_matmul(a, param)
              lhs <- -1/2 * a + torch::torch_eye(nrow(a), ncol(a),
                                                 dtype = torch::torch_double())

              new_param <- torch::linalg_solve(lhs, rhs)

              if (!is.na(group$threshold)){
                new_param <- soft(new_param, group$threshold)
                res_svd <- torch::linalg_svd(new_param, full_matrices=FALSE)
                proj_param <- torch::torch_matmul(res_svd[[1]], res_svd[[3]])
                param$copy_(proj_param)
              } else {
                param$copy_(new_param)
              }






            }
          }
       else {
          for (p in seq_along(group$params)) {
            param <- group$params[[p]]

            if (is.null(param$grad) || torch::is_undefined_tensor(param$grad)) {
              next
            }

            param$add_(param$grad, alpha = -group$lr)
          }
        }
      }

    })
  }
)
