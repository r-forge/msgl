\name{sgl_predict}
\alias{sgl_predict}
\title{Sgl predict}
\usage{
  sgl_predict(module_name, PACKAGE, object, data, ...)
}
\arguments{
  \item{module_name}{reference to objective specific C++
  routines}

  \item{PACKAGE}{name of the package}

  \item{data}{a list of data objects -- will be parsed to
  the specified module}

  \item{...}{}
}
\value{
  sgl object content will depend on the C++ response class
}
\description{
  Sgl predict
}
\author{
  Martin Vincent
}

