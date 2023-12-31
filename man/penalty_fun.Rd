% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{penalty_fun}
\alias{penalty_fun}
\title{Penalty constraints for robust and sparse PCA / CCA}
\usage{
penalty_fun(z, pen, alpha)
}
\arguments{
\item{z}{A torch::torch_tensor}

\item{pen}{penalty parameter}

\item{alpha}{elastic net mixing parameter}
}
\value{
A 1x1 torch::torch_tensor
}
\description{
Function to compute elastic net penalty for a tensor
}
