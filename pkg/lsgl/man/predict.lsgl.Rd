\name{predict.lsgl}
\alias{predict.lsgl}
\title{Predict}
\usage{
  predict.lsgl(object, x, ...)
}
\arguments{
  \item{object}{}

  \item{x}{}

  \item{...}{}
}
\value{
  \item{Yhat}{the linear predictors -- a list of length
  \code{length(fit$beta)} one item for each model, with
  each item a matrix of size \eqn{N_\textrm{new} \times K}
  containing the linear predictors.}
}
\description{
  Predict
}
\author{
  Martin Vincent
}

