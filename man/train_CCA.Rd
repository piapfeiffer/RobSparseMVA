% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_CCA_MM.R
\name{train_CCA}
\alias{train_CCA}
\title{MM algorithm for CCA objects
implements MM algorithms for robust and sparse CCA}
\usage{
train_CCA(
  C,
  p,
  q,
  pen_x,
  pen_y,
  order,
  alpha_x,
  alpha_y,
  tol,
  lr,
  epochs,
  low_a,
  low_b,
  lr_decay,
  show_warnings = FALSE,
  init_ortho = TRUE
)
}
\arguments{
\item{C}{estimated covariance matrix}

\item{p}{number of variables in data_x}

\item{q}{number of variables in data_y}

\item{pen_x}{(1xk) A positive real vector containing the penalty parameter for
data_x for association of order k.}

\item{pen_y}{(1xk) A positive real vector containing the penalty parameter for
data_y for association of order k.}

\item{order}{The order of association to be computed}

\item{alpha_x}{(1xk) A positive real vector containing numbers between 0 and
1, indicating the elastic
 net parameter for data_x for association of order k.}

\item{alpha_y}{(1xk) A positive real vector containing numbers between 0 and
1, indicating the elastic
 net parameter for data_y for association of order k.}

\item{tol}{tolerance for training algorithm}

\item{lr}{learning rate for gradient descent training algorithm}

\item{epochs}{maximum number of epochs for gradient descent training algorithm}

\item{low_a}{A (px(k-1)) matrix of lower order linear combinations}

\item{low_b}{A (qx(k-1)) matrix of lower order linear combinations}

\item{lr_decay}{parameter for}

\item{show_warnings}{logical}

\item{init_ortho}{logical}
}
\description{
MM algorithm for CCA objects
implements MM algorithms for robust and sparse CCA
}
