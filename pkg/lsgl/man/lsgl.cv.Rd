\name{lsgl.cv}
\alias{lsgl.cv}
\title{Cross validation using multiple possessors}
\usage{
  lsgl.cv(x, y, intercept = TRUE,
    covariateGrouping = factor(1:ncol(x)),
    groupWeights = c(sqrt(ncol(y) * table(covariateGrouping))),
    parameterWeights = matrix(1, nrow = ncol(y), ncol = ncol(x)),
    alpha = 0.5, lambda, fold = 10L, cv.indices = list(),
    max.threads = 2L,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{}

  \item{y}{}

  \item{intercept}{}

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

