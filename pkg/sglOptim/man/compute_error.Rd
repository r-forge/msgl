\name{compute_error}
\alias{compute_error}
\title{Helper function for computing error rates}
\usage{
  compute_error(object, data = NULL, response.name,
    response, loss)
}
\arguments{
  \item{object}{a object containing responses}

  \item{data}{a data object}

  \item{response.name}{the name of the response, if
  \code{response.name = NULL} then \code{x} will be treated
  as the response.}

  \item{response}{the response}

  \item{loss}{the loss function}
}
\value{
  a vector with the computed error rates
}
\description{
  Helper function for computing error rates
}
\details{
  This function can be used to compute error rates. It is
  consist with the use cases of the \code{Err} genetic
  function. (see \pkg{msgl} package for an example of how
  to use this function)
}
\author{
  Martin Vincent
}

