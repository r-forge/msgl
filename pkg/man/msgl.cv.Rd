\name{msgl.cv}
\alias{msgl.cv}
\title{Multinomial sparse group lasso cross validation}
\usage{
  msgl.cv(x, classes,
    featureWeights = .featureWeights(x, classes),
    classWeights = .classWeights(classes), alpha = 0.5,
    standardize = TRUE, lambda, fold = 10L,
    cv.classes = .list.factor(classes),
    cv.indices = list(), do.refit = FALSE,
    sparse.data = FALSE, max.threads = 2L, seed = 331L,
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

  \item{standardize}{if TRUE the covariates are standardize
  before fitting the model. The model parameters are
  returned in the original scale.}

  \item{lambda}{the lambda sequence for the regularization
  path}

  \item{fold}{}

  \item{cv.classes}{}

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
  \item{features}{Average number of features used in the
  models} \item{parameters}{Average number of parameters
  used in the models}
}
\description{
  Multinomial sparse group lasso cross validation
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)

# Missclassification count
colSums(fit.cv$classes != classes)
}
\author{
  Martin Vincent
}

