\name{predict.lrsgl}
\alias{predict.lrsgl}
\title{Predict}
\usage{
  predict.lrsgl(fit, x, sparse.data = FALSE)
}
\arguments{
  \item{fit}{object of class 'lrsgl'}

  \item{x}{new data}

  \item{sparse.data}{if TRUE x will be treated as sparse}
}
\value{
  \item{link}{linear predictors} \item{response}{estimated
  probabilities} \item{classes}{estimated classes}
}
\description{
  Predict
}
\author{
  Martin Vincent
}

