\name{sgl_fit}
\alias{sgl_fit}
\title{Fit a sparse group lasso regularization path.}
\usage{
  sgl_fit(module_name, PACKAGE, data, covariateGrouping,
    groupWeights, parameterWeights, alpha, lambda,
    return = 1:length(lambda),
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{call_sym}{reference to objective specific C++
  routines}

  \item{data}{}

  \item{covariateGrouping}{grouping of covariates, a vector
  of length \eqn{p}. Each element of the vector specifying
  the group of the covariate.}

  \item{groupWeights}{the group weights, a vector of length
  \eqn{m} (the number of groups).}

  \item{parameterWeights}{a matrix of size \eqn{K \times
  (p)}.}

  \item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1
  for lasso, between 0 and 1 gives a sparse group lasso
  penalty.}

  \item{lambda}{the lambda sequence for the regularization
  path.}

  \item{return}{the indices of lambda values for which to
  return a the fitted parameters.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  \item{beta}{the fitted parameters -- a list of length
  \code{length(lambda)} with each entry a matrix of size
  \eqn{K\times (p+1)} holding the fitted parameters}
  \item{loss}{the values of the loss function}
  \item{objective}{the values of the objective function
  (i.e. loss + penalty)} \item{lambda}{the lambda values
  used}
}
\description{
  Fit a sparse group lasso regularization path.
}
\author{
  Martin Vincent
}

