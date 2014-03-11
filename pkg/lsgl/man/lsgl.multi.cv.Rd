\name{lsgl.multi.cv}
\alias{lsgl.multi.cv}
\title{Cross validation using multiple possessors}
\usage{
  lsgl.multi.cv(X, Y, weights = rep(1/nrow(X), nrow(X)),
    covariateGrouping = factor(1:ncol(X)),
    groupWeights = NULL, parameterWeights = NULL,
    alpha = 0.5, lambda, fold = 10L, cv.indices = list(),
    max.threads = 2L,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{X}{design matrix, matrix of size \eqn{N \times p}.}

  \item{Y}{response matrix, matrix of size \eqn{N \times
  K}.}

  \item{weights,}{sample weights a vector of size \eqn{N}.}

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

  \item{lambda}{lambda sequence.}

  \item{fold}{}

  \item{cv.indices}{}

  \item{max.threads}{}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  sgl object.
}
\description{
  Cross validation using multiple possessors
}
\author{
  Martin Vincent
}

