\name{sgl_predict}
\alias{sgl_predict}
\title{Sgl predict}
\usage{
  sgl_predict(module_name, PACKAGE, object, data, ...)
}
\arguments{
  \item{module_name}{reference to objective specific C++
  routines.}

  \item{PACKAGE}{name of the calling package.}

  \item{object}{a sgl object containing a list of estimated
  models.}

  \item{data}{a list of data objects -- will be parsed to
  the specified module.}

  \item{...}{not used.}
}
\value{
  \item{responses}{content will depend on the C++ response
  class} \item{lambda}{the lambda sequence used.}
}
\description{
  Sgl predict
}
\author{
  Martin Vincent
}

