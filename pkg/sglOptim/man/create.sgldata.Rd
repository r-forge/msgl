\name{create.sgldata}
\alias{create.sgldata}
\title{Create sgl data}
\usage{
  create.sgldata(x, y, weights = rep(1/nrow(x), nrow(x)),
    sampleGrouping = NULL)
}
\arguments{
  \item{x}{design matrix, matrix of size \eqn{N \times p}.}

  \item{y}{responses, vector of length \eqn{N}.}

  \item{sampleWeights}{sample weights, a vector of length
  \eqn{N}.}

  \item{sampleGrouping}{grouping of samples, a factor of
  length \eqn{N}. Default is no grouping (NULL), that is
  all samples is the same group}
}
\description{
  Create sgl data
}
\author{
  Martin Vincent
}

