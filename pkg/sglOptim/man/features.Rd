\name{features}
\alias{features}
\title{Generic function for extracting nonzero features (or groups)}
\usage{
  features(object, ...)
}
\arguments{
  \item{object}{an object}

  \item{...}{additional parameters (optional)}
}
\value{
  a list of length \code{nmod(x)} containing the nonzero
  features of the models.
}
\description{
  Extracts nonzero features for each model.
}
\author{
  Martin Vincent
}

