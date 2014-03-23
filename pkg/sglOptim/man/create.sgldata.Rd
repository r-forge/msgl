\name{create.sgldata}
\alias{create.sgldata}
\title{Create a sgldata object}
\usage{
  create.sgldata(x, y, weights = rep(1/nrow(x), nrow(x)),
    sampleGrouping = NULL, group.names = NULL,
    sparseX = is(x, "sparseMatrix"))
}
\arguments{
  \item{x}{the design matrix, a matrix of size \eqn{N
  \times p} (will be parsed to the loss module as
  \code{X}).}

  \item{y}{the responses, \code{NULL}, a vector or a matrix
  (will be parsed to the loss module as \code{Y})..}

  \item{weights}{sample weights, a vector of length \eqn{N}
  (will be parsed to the loss module as \code{W}).}

  \item{sampleGrouping}{grouping of samples, a factor of
  length \eqn{N} (will be parsed to the loss module as
  \code{G}). Default is no grouping (NULL), that is all
  samples is the same group.}

  \item{group.names}{a vector with the names of the
  parameter groups (the length must equal the number of
  rows in the \eqn{\beta} matrix).}

  \item{sparseX}{if TRUE \code{x} will be treated as
  sparse, if FALSE \code{x} will be treated as dens.}
}
\description{
  Creates a sgldata object from a design matrix and an
  optional response vector or matrix.
}
\author{
  Martin Vincent
}
\seealso{
  Other sgldata: \code{\link{prepare.args}},
  \code{\link{prepare.args.sgldata}},
  \code{\link{rearrange.sgldata}}
}

