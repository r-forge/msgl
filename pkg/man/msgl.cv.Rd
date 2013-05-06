\name{msgl.cv}
\alias{msgl.cv}
\title{Multinomial sparse group lasso cross validation using multiple possessors}
\usage{
  msgl.cv(x, classes, sampleWeights = NULL,
    grouping = NULL, groupWeights = NULL,
    parameterWeights = NULL, alpha = 0.5,
    standardize = TRUE, lambda, fold = 10L,
    cv.indices = list(), sparse.data = FALSE,
    max.threads = 2L, seed = 331L,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{Design matrix, matrix of size \eqn{N \times p}.}

  \item{classes}{Classes, factor of length \eqn{N}.}

  \item{sampleWeights}{Sample weights, a vector of length
  \eqn{N}.}

  \item{grouping}{Grouping of covariates, a vector of
  length \eqn{p}. Each element of the vector specifying the
  group of the covariate.}

  \item{groupWeights}{The group weights, a vector of length
  \eqn{m+1} (the number of groups). The first element of
  the vector is the intercept weight. If \code{groupWeights
  = NULL} default weights will be used. Default weights are
  0 for the intercept and \eqn{\sqrt{K\cdot\textrm{number
  of covariates in the group}}} for all other weights.}

  \item{parameterWeights}{A matrix of size \eqn{K \times
  (p+1)}. The first column of the matrix is the intercept
  weights. Default weights are is 0 for the intercept
  weights and 1 for all other weights.}

  \item{alpha}{Alpha value 0 for group lasso, 1 for lasso,
  between 0 and 1 gives a sparse group lasso penalty.}

  \item{standardize}{If TRUE the covariates are standardize
  before fitting the model. The model parameters are
  returned in the original scale.}

  \item{lambda}{The lambda sequence for the regularization
  path.}

  \item{fold}{the fold of the cross validation. The data
  will be split into fold disjoint subsets, keeping the
  ration of classes approximately equal.}

  \item{cv.indices}{a list of indices of a cross validation
  splitting. If \code{cv.indices = NULL} then a random
  splitting will be generated using the \code{fold}
  argument.}

  \item{sparse.data}{If TRUE \code{x} will be treated as
  sparse, if \code{x} is a sparse matrix it will be treated
  as sparse by default.}

  \item{max.threads}{the maximal number of threads to be
  used}

  \item{seed}{the seed used for generating the random cross
  validation splitting, only used if \code{cv.indices =
  NULL}.}

  \item{algorithm.config}{The algorithm configuration to be
  used.}
}
\value{
  \item{link}{The linear predictors - a list of length
  \code{length(lambda)} one item for each lambda value,
  with each item a matrix of size \eqn{K \times N}
  containing the linear predictors} \item{response}{The
  estimated probabilities - a list of length
  \code{length(lambda)} one item for each lambda value,
  with each item a matrix of size \eqn{K \times N}
  containing the probabilities} \item{classes}{The
  estimated classes - a matrix of size \eqn{N \times d}
  with \eqn{d=}\code{length(lambda)}} \item{cv.indices}{The
  cross validation splitting used} \item{features}{Average
  number of features used in the models}
  \item{parameters}{Average number of parameters used in
  the models}
}
\description{
  Multinomial sparse group lasso cross validation using
  multiple possessors
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

