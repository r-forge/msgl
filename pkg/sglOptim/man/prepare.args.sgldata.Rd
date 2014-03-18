\name{prepare.args.sgldata}
\alias{prepare.args.sgldata}
\title{Prepare sgl function arguments}
\usage{
  \method{prepare.args}{sgldata} (data, parameterGrouping,
    groupWeights, parameterWeights, alpha, ...)
}
\arguments{
  \item{data}{sgldata object}

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

  \item{...}{not used}
}
\description{
  Prepare sgl function arguments
}
\author{
  Martin Vincent
}
\seealso{
  Other sgldata: \code{\link{create.sgldata}},
  \code{\link{rearrange.sgldata}}
}

