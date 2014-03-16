\name{coef}
\alias{coef}
\title{Generic function for extracting the nonzero coefficients}
\usage{
  coef(x, index, ...)
}
\arguments{
  \item{x}{an object}

  \item{index}{a vector of indices of the models}

  \item{...}{additional paramters (optional)}
}
\value{
  a list of length \code{length(index)} containing the
  nonzero paramters of the models
}
\description{
  Generic function for extracting the nonzero coefficients
}
\author{
  Martin Vincent
}

