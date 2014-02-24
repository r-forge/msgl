\name{sgl_lambda_sequence}
\alias{sgl_lambda_sequence}
\title{Generic routine for computing a lambda sequence for the regularization path}
\usage{
  sgl_lambda_sequence(module_name, PACKAGE, data,
    parameterGrouping, groupWeights, parameterWeights,
    alpha = 0.5, d = 100L, lambda.min,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{module_name}{reference to objective specific C++
  routines.}

  \item{PACKAGE}{name of the calling package.}

  \item{data}{a list of data objects -- will be parsed to
  the specified module.}

  \item{parameterGrouping}{grouping of parameters, a vector
  of length \eqn{p}. Each element of the vector specifying
  the group of the parameters in the corresponding column
  of \eqn{\beta}.}

  \item{groupWeights}{the group weights, a vector of length
  \code{length(unique(parameterGrouping))} (the number of
  groups).}

  \item{parameterWeights}{a matrix of size \eqn{q \times
  p}.}

  \item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1
  for lasso, between 0 and 1 gives a sparse group lasso
  penalty.}

  \item{d}{the length of lambda sequence.}

  \item{lambda.min}{the smallest lambda value in the
  computed sequence.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  a vector of length \code{d} containing the compute lambda
  sequence.
}
\description{
  Computes a decreasing lambda sequence of length \code{d}.
  The sequence ranges from a data determined maximal lambda
  \eqn{\lambda_\textrm{max}} to the user inputed
  \code{lambda.min}.
}
\author{
  Martin Vincent
}

