\name{sgl_cv}
\alias{sgl_cv}
\title{Sparse group lasso cross validation using multiple possessors}
\usage{
  sgl_cv(module_name, PACKAGE, data, covariateGrouping,
    groupWeights, parameterWeights, alpha, lambda,
    fold = 2L, cv.indices = list(), max.threads = 2L,
    seed = 331L, algorithm.config = sgl.standard.config)
}
\arguments{
  \item{call_sym}{reference to objective specific C++
  routines}

  \item{data}{}

  \item{covariateGrouping}{grouping of covariates, a vector
  of length \eqn{p}. Each element of the vector specifying
  the group of the covariate.}

  \item{groupWeights}{the group weights, a vector of length
  \eqn{m+1} (the number of groups).}

  \item{parameterWeights}{a matrix of size \eqn{K \times
  (p+1)}.}

  \item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1
  for lasso, between 0 and 1 gives a sparse group lasso
  penalty.}

  \item{lambda}{the lambda sequence for the regularization
  path.}

  \item{fold}{the fold of the cross validation, an integer
  larger than \eqn{1} and less than \eqn{N+1}. Ignored if
  \code{cv.indices != NULL}. If
  \code{fold}\eqn{\le}\code{max(table(classes))} then the
  data will be split into \code{fold} disjoint subsets
  keeping the ration of classes approximately equal.
  Otherwise the data will be split into \code{fold}
  disjoint subsets without keeping the ration fixed.}

  \item{cv.indices}{a list of indices of a cross validation
  splitting. If \code{cv.indices = NULL} then a random
  splitting will be generated using the \code{fold}
  argument.}

  \item{max.threads}{the maximal number of threads to be
  used}

  \item{seed}{the seed used for generating the random cross
  validation splitting, only used if
  \code{fold}\eqn{\le}\code{max(table(classes))}.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\description{
  Sparse group lasso cross validation using multiple
  possessors
}
\author{
  Martin Vincent
}

