\name{lsgl}
\alias{lsgl}
\title{Fit}
\usage{
  lsgl(x, y, weights = rep(1/nrow(x), nrow(x)),
    sampleGrouping = factor(rep(1, nrow(x))),
    covariateGrouping = factor(1:ncol(x)),
    groupWeights = c(sqrt(length(levels(sampleGrouping)) * table(covariateGrouping))),
    parameterWeights = matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x)),
    alpha = 0.5, lambda,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{}

  \item{y}{}

  \item{weights}{}

  \item{sampleGrouping}{}

  \item{covariateGrouping}{grouping of covariates, a factor
  or vector of length \eqn{p}. Each element of the
  factor/vector specifying the group of the covariate.}

  \item{groupWeights}{the group weights, a vector of length
  \eqn{m+1} (the number of groups).}

  \item{parameterWeights}{a matrix of size \eqn{K \times
  (p+1)}.}

  \item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1
  for lasso, between 0 and 1 gives a sparse group lasso
  penalty.}

  \item{lambd}{lambda sequence.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  a vector of length \code{d} containing the compute lambda
  sequence.
}
\description{
  Fit
}
\author{
  Martin Vincent
}

