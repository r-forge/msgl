\name{predict.msgl}
\alias{predict.msgl}
\title{Predict}
\usage{
  predict.msgl(object, x, sparse.data = FALSE, ...)
}
\arguments{
  \item{object}{an object of class msgl, produced with
  \code{msgl}.}

  \item{x}{a data matrix of size \eqn{N_\textrm{new} \times
  p}.}

  \item{sparse.data}{If TRUE \code{x} will be treated as
  sparse, if \code{x} is a sparse matrix it will be treated
  as sparse by default.}
}
\value{
  \item{link}{The linear predictors - a list of length
  \code{length(fit$beta)} one item for each model, with
  each item a matrix of size \eqn{K \times N_\textrm{new}}
  containing the linear predictors} \item{response}{The
  estimated probabilities - a list of length
  \code{length(fit$beta)} one item for each model, with
  each item a matrix of size \eqn{K \times N_\textrm{new}}
  containing the probabilities} \item{classes}{The
  estimated classes - a matrix of size \eqn{N_\textrm{new}
  \times d} with \eqn{d=}\code{length(fit$beta)}}
}
\description{
  Computes the linear predictors, the estimated
  probabilities and the estimated classes for a new data
  set.
}
\author{
  Martin Vincent
}

