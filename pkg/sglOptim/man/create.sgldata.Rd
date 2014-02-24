\name{create.sgldata}
\alias{create.sgldata}
\title{Create sgl data}
\usage{
  create.sgldata(x, y, weights = rep(1/nrow(x), nrow(x)),
    sampleGrouping = NULL, sparseX = is(x, "sparseMatrix"))
}
\arguments{
  \item{x}{design matrix, matrix of size \eqn{N \times p}.}

  \item{y}{responses, vector of length \eqn{N}.}

  \item{weights}{sample weights, a vector of length
  \eqn{N}.}

  \item{sampleGrouping}{grouping of samples, a factor of
  length \eqn{N}. Default is no grouping (NULL), that is
  all samples is the same group.}

  \item{sparseX}{if TRUE \code{x} will be treated as
  sparse, if FALSE \code{x} will be treated as dens.}
}
\description{
  Create sgl data
}
\author{
  Martin Vincent
}

