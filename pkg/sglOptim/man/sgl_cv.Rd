\name{sgl_cv}
\alias{sgl_cv}
\title{Generic sparse group lasso cross validation using multiple possessors}
\usage{
  sgl_cv(module_name, PACKAGE, data, parameterGrouping,
    groupWeights, parameterWeights, alpha, lambda,
    fold = 2, cv.indices = list(), max.threads = 2,
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
  used.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  \item{responses}{content will depend on the C++ response
  class} \item{cv.indices}{the cross validation splitting
  used} \item{features}{number of features used in the
  models} \item{parameters}{number of parameters used in
  the models} \item{lambda}{the lambda sequence used.}
}
\description{
  Generic sparse group lasso cross validation using
  multiple possessors
}
\author{
  Martin Vincent
}

