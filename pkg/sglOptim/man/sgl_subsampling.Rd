\name{sgl_subsampling}
\alias{sgl_subsampling}
\title{Generic sparse group lasso subsampling procedure}
\usage{
  sgl_subsampling(module_name, PACKAGE, data,
    parameterGrouping, groupWeights, parameterWeights,
    alpha, lambda, training, test, max.threads = 2L,
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

  \item{training}{a list of training samples, each item of
  the list corresponding to a subsample. Each item in the
  list must be a vector with the indices of the training
  samples for the corresponding subsample. The length of
  the list must equal the length of the \code{test} list.}

  \item{test}{a list of test samples, each item of the list
  corresponding to a subsample. Each item in the list must
  be vector with the indices of the test samples for the
  corresponding subsample. The length of the list must
  equal the length of the \code{training} list.}

  \item{max.threads}{the maximal number of threads to be
  used.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  sgl object content will depend on the C++ response class
}
\description{
  Support the use of multiple processors.
}
\author{
  Martin Vincent
}

