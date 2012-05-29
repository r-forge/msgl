\name{msgl.subsampleing}
\alias{msgl.subsampleing}
\title{Multinomial sparse group lasso subsampling}
\usage{
  msgl.subsampleing(x, classes,
    featureWeights = .featureWeights(x, classes),
    classWeights = .classWeights(classes), alpha = 0.5,
    standardize = TRUE, lambda,
    cv.classes = .list.factor(classes), d = 100L,
    fraction = 0.2, subsamples = list(), do.refit = FALSE,
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

  \item{cv.classes}{unused}

  \item{d}{number of subsamples}

  \item{fraction}{}

  \item{subsamples}{alternative way of specifing the
  subsamples. A list of index vectors.}

  \item{do.refit}{if TRUE a refitted model will be
  returned}

  \item{sparse.data}{if TRUE x will be treated as sparse}

  \item{max.threads}{maximal number of threads}

  \item{seed}{}

  \item{algorithm.config}{the algorithm configuration to be
  used}
}
\value{
  ...
}
\description{
  Multinomial sparse group lasso subsampling
}
\author{
  Martin vincent
}

