# Hyperparameter Optimization for routines in RobSparseMVA

#' Bayesian Hyperparameter optimization
#' Get optimal penalties using bayesian optimization.
#' Score function depends on whether sparsity is assumed or not
#' @param C estimated covariance matrix
#' @param p number of variables in data_x
#' @param q number of variables in data_y
#' @param n number of observations
#' @param alpha_x (1xk) A positive real vector containing numbers between 0 and
#' 1, indicating the elastic
#'  net parameter for data_x for association of order k.
#' @param alpha_y (1xk) A positive real vector containing numbers between 0 and
#' 1, indicating the elastic
#'  net parameter for data_y for association of order k.
#' @param order The order of association to be computed
#' @param low_a A (px(k-1)) matrix of lower order linear combinations
#' @param low_b A (qx(k-1)) matrix of lower order linear combinations
#' @param tol tolerance for training algorithm
#' @param lr  learning rate for gradient descent training algorithm
#' @param epochs maximum number of epochs for gradient descent training algorithm
#' @param lr_decay parameter for
#' @param criterion type of criterion used to determine optimal hyperparameters,
#' default: criterion = TPO
#' @param init_ortho logical
#' @returns a list object containing the optimal hyperparameters and a summary
#' of the optimization process
bayesian_optimization_CCA <- function(C, p, q, n,
                                      alpha_x, alpha_y,
                                      order,
                                      low_a, low_b,
                                      tol, lr, epochs,
                                      lr_decay,
                                      criterion = "TPO",
                                      init_ortho = TRUE){

  scoringFunction <- function(pen_x, pen_y){
    res_train <- train_CCA(C, p, q,
                           pen_x * sqrt(min(p, n)),
                           pen_y * sqrt(min(q, n)),
                           order,
                           alpha_x[order], alpha_y[order],
                           tol = tol, lr = lr, epochs = epochs,
                           low_a = low_a,
                           low_b = low_b,
                           lr_decay = lr_decay,
                           init_ortho = init_ortho)

    R <- res_train$best_measure
    L <- res_train$best_loss
    l1_a <- norm(res_train$best_a, "1")
    l1_b <- norm(res_train$best_b, "1")
    p_rate <- sum(res_train$best_a!=0)
    q_rate <- sum(res_train$best_b!=0)

    if (criterion == "TPO"){
      score <- abs(R)*(2 - p_rate/p * alpha_x[order] - q_rate/q * alpha_y[order] + 1e-5)
    } else if (criterion == "BIC"){
      rlang::abort('Not yet implemented', class = "feature_error")
    } else {
      score <- abs(R)
    }

    pred <- abs(R)

    return(list(Score = score, Pred = pred,
                Sparse_x = 1-p_rate/p,
                Sparse_y = 1-q_rate/q
                ))
  }

  bounds <- list(pen_x = c(1 / sqrt(min(p, n)), 1),
                 pen_y = c(1 / sqrt(min(q, n)), 1))

  initGrid <- data.frame(pen_x = seq(from = min(bounds$pen_x) ,to = max(bounds$pen_x), length.out = 5),
                         pen_y = seq(from = min(bounds$pen_y) ,to = max(bounds$pen_y), length.out = 5))

  summary <- c()
  optObj <- c()
  optObj$scoreSummary <- NA
  best_params <- list(pen_x = max(bounds$pen_x)* sqrt(min(p, n)), pen_y = max(bounds$pen_y)* sqrt(min(q, n)))
  summary_list <- list()

  out <- tryCatch({

    optObj <- ParBayesianOptimization::bayesOpt(
      FUN = scoringFunction
      , bounds = bounds
      , acq = "ei"
      , initGrid = initGrid
      , iters.n = 10
      , eps = 0.1
      , parallel = FALSE
      , errorHandling = 3
      , acqThresh = 1
      , plotProgress = FALSE
      , verbose = 1
    )
    optObj$scoreSummary$pen_x <- optObj$scoreSummary$pen_x * sqrt(min(p, n))
    optObj$scoreSummary$pen_y <- optObj$scoreSummary$pen_y * sqrt(min(q, n))

    best_params <- ParBayesianOptimization::getBestPars(optObj)
  },  error = function(err) {
    print(paste("ERROR:  ",err, "\n maximum penalty values are used!"))
  }, finally = {
    return(list(best_params = best_params, summary = optObj$scoreSummary))
  })
}
