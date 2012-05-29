\name{msgl.lambda.seq}
\alias{msgl.lambda.seq}
\title{Computes a lambda sequence for the regularization path}
\usage{
  msgl.lambda.seq(x, classes,
    featureWeights = .featureWeights(x, classes),
    classWeights = .classWeights(classes), alpha = 0.5,
    d = 100L, standardize = TRUE, lambda.min,
    sparse.data = FALSE,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{design matrix, matrix of size N x p}

  \item{classes}{grouping of samples, factor of length N}

  \item{featureWeights}{the group weights, a vector of
  length p}

  \item{classWeights}{a vector of length K}

  \item{alpha}{0 for group lasso, 1 for lasso, between 0
  and 1 gives a sparse group lasso penalty}

  \item{d}{length of the required sequence}

  \item{standardize}{if TRUE the covariates are standardize
  before fitting the model. The model parameters are
  returned in the original scale.}

  \item{lambda.min}{minimum lambda value}

  \item{sparse.data}{if TRUE x will be treated as sparse}

  \item{algorithm.config}{the algorithm configuration to be
  used}
}
\value{
  lambda sequence
}
\description{
  Computes a lambda sequence for the regularization path
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
}
\author{
  Martin Vincent
}

