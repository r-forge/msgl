\name{logreg.sgl}
\alias{logreg.sgl}
\title{Fit a logistic sparse group lasso regularization path}
\usage{
  logreg.sgl(x, classes, group.dim, groupWeights,
    parameterWeights, alpha = 0.5, standardize = TRUE,
    lambda, return = 1:length(lambda), do.refit = FALSE,
    sparse.data = FALSE,
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

  \item{standardize}{if TRUE the covariates are standardize
  before fitting the model. The model parameters are
  returned in the original scale.}

  \item{lambda}{the lambda sequence for the regularization
  path}

  \item{return}{the indices of lambda values for which to
  return a the fitted parameters}

  \item{do.refit}{if TRUE a refitted model will be
  returned}

  \item{sparse.data}{if TRUE x will be treated as sparse}

  \item{algorithm.config}{the algorithm configuration to be
  used}
}
\value{
  \item{beta}{the fittede paramters} \item{loss}{the values
  of the loss function} \item{objective}{the values of the
  objective function (i.e. loss + penalty)}
  \item{lambda}{the lambda values used}
  \item{beta.refit}{refitted paramters (only if do_refit =
  TRUE)} \item{loss.refit}{the values of the refitted loss
  function (only if do_refit = TRUE)}
}
\description{
  Note that, if standardize = TRUE the fitted parameters
  are currently not rescaled to orgnial scale before
  returned.
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
fit <- logreg.sgl(x, classes, group.dim, group.weights, parameter.weights, lambda = lambda, standardize = FALSE)
fit$beta[[10]] #model with lambda = lambda[10]
}
\author{
  Martin Vincent
}

