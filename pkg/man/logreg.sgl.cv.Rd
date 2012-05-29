\name{logreg.sgl.cv}
\alias{logreg.sgl.cv}
\title{Logistic sparse group lasso cross validation}
\usage{
  logreg.sgl.cv(x, classes, group.dim, groupWeights,
    parameterWeights, alpha = 0.5, standardize = TRUE,
    lambda, fold = 10L, cv.indices = list(),
    do.refit = FALSE, sparse.data = FALSE,
    max.threads = 2L, seed = 331L,
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

  \item{fold}{}

  \item{cv.indices}{}

  \item{do.refit}{if TRUE a refitted model will be
  returned}

  \item{sparse.data}{if TRUE x will be treated as sparse}

  \item{max.threads}{maximal number of threads}

  \item{seed}{}

  \item{algorithm.config}{the algorithm configuration to be
  used}
}
\value{
  \item{link}{linear predictors} \item{response}{estimated
  probabilities} \item{classes}{estimated classes}
  \item{link.refit}{linear predictors for the refitted
  models (only if do_refit = TRUE)}
  \item{response.refit}{estimated probabilities for the
  refitted models (only if do_refit = TRUE)}
  \item{classes.refit}{estimated classes for the refitted
  models (only if do_refit = TRUE)} \item{cv.indices}{}
  \item{groups}{Average number of groups used in the
  models} \item{parameters}{Average number of parameters
  used in the models}
}
\description{
  Logistic sparse group lasso cross validation
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
fit.cv <- logreg.sgl.cv(x, classes, group.dim, group.weights, parameter.weights, lambda = lambda, fold = 10L, max.threads = 2L)

# Missclassification count
colSums(fit.cv$classes != classes)
}
\author{
  Martin Vincent
}

