\name{models.sgl}
\alias{models.sgl}
\title{Returns the estimated models (that is the \eqn{beta} matrices)}
\usage{
  \method{models}{sgl} (object, index = 1:nmod(object),
    ...)
}
\arguments{
  \item{object}{a sgl object}

  \item{index}{indices of the models to be returned}

  \item{...}{not used}
}
\value{
  a list of sparse matrices
}
\description{
  Returns the estimated models (that is the \eqn{beta}
  matrices)
}
\author{
  Martin Vincent
}

