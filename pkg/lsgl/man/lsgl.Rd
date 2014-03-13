\name{lsgl}
\alias{lsgl}
\title{Fit a linear multi-response model using sparse group lasso}
\usage{
  lsgl(x, y, intercept = TRUE,
    covariateGrouping = factor(1:ncol(x)),
    groupWeights = c(sqrt(ncol(y) * table(covariateGrouping))),
    parameterWeights = matrix(1, nrow = ncol(y), ncol = ncol(x)),
    alpha = 0.5, lambda,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{design matrix, matrix of size \eqn{N \times p}.}

  \item{y}{response matrix, matrix of size \eqn{N \times
  K}.}

  \item{intercept}{}

  \item{covariateGrouping}{grouping of covariates, a factor
  or vector of length \eqn{p}. Each element of the
  factor/vector specifying the group of the covariate.}

  \item{groupWeights}{the group weights, a vector of length
  \eqn{m+1} (the number of groups).}

  \item{parameterWeights}{a matrix of size \eqn{K \times
  (p+1)}.}

  \item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1
  for lasso, between 0 and 1 gives a sparse group lasso
  penalty.}

  \item{lambda}{lambda sequence.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  \item{beta}{the fitted parameters -- a list of length
  \code{length(return)} with each entry a matrix of size
  \eqn{q\times (p+1)} holding the fitted parameters.}
  \item{loss}{the values of the loss function.}
  \item{objective}{the values of the objective function
  (i.e. loss + penalty).} \item{lambda}{the lambda values
  used.}
}
\description{
  For a linear multi-response model with \eqn{p} covariates
  dived into \eqn{m} groups using sparse group lasso.
}
\details{
  This function computes a sequence of minimizers (one for
  each lambda given in the \code{lambda} argument) of
  \deqn{\|Y-X\beta\|_F^2 + \lambda \left( (1-\alpha)
  \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha
  \sum_{i=1}^{n} \xi_i |\beta_i| \right)} where
  \eqn{\|\cdot\|_F} is the frobenius norm. The vector
  \eqn{\beta^{(J)}} denotes the parameters associated with
  the \eqn{J}'th group of covariates. The group weights are
  denoted by \eqn{\gamma \in [0,\infty)^m} and the
  parameter weights by \eqn{\xi = (\xi^{(1)},\dots,
  \xi^{(m)}) \in [0,\infty)^n} with \eqn{\xi^{(1)}\in
  [0,\infty)^{n_1},\dots, \xi^{(m)} \in [0,\infty)^{n_m}}.
}
\author{
  Martin Vincent
}

