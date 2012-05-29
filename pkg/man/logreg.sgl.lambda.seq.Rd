\name{logreg.sgl.lambda.seq}
\alias{logreg.sgl.lambda.seq}
\title{Computes a lambda sequence for the logistic sparse group lasso regularization path}
\usage{
  logreg.sgl.lambda.seq(x, classes, group.dim,
    groupWeights, parameterWeights, alpha = 0.5, d = 100L,
    standardize = TRUE, lambda.min, sparse.data = FALSE,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{design matrix, matrix of size N x p}

  \item{classes}{grouping of samples, factor of length N}

  \item{group.dim}{vector of group dimensions}

  \item{groupWeights}{group weights}

  \item{parameterWeights}{parameter weights}

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
  Computes a lambda sequence for the logistic sparse group
  lasso regularization path
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes > 5

group.dim <- as.integer(c(1, rep(10, 40)))
group.dim <- as.integer(group.dim)
group.weights <- rep(1, length(group.dim))
parameter.weights <- rep(1, 401)
group.weights[1] <- 0 # no penalty on the intercept
parameter.weights[1] <- 0

lambda <- logreg.sgl.lambda.seq(x, classes, group.dim, group.weights, parameter.weights, lambda.min = 0.01)
}
\author{
  Martin Vincent
}

