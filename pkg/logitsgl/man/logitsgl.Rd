% Generated by roxygen2 (4.0.1): do not edit by hand
\name{logitsgl}
\alias{logitsgl}
\title{Fit a linear multiple output model using sparse group lasso}
\usage{
logitsgl(x, y, intercept = TRUE, grouping = factor(1:ncol(x)),
  groupWeights = c(sqrt(ncol(y) * table(grouping))),
  parameterWeights = matrix(1, nrow = ncol(y), ncol = ncol(x)), alpha = 1,
  lambda, algorithm.config = logitsgl.standard.config)
}
\arguments{
\item{x}{design matrix, matrix of size \eqn{N \times p}.}

\item{y}{response matrix, matrix of size \eqn{N \times K}.}

\item{intercept}{should the model include intercept parameters}

\item{grouping}{grouping of features, a factor or vector of length \eqn{p}.
Each element of the factor/vector specifying the group of the feature.}

\item{groupWeights}{the group weights, a vector of length \eqn{m} (the number of groups).}

\item{parameterWeights}{a matrix of size \eqn{K \times p}.}

\item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.}

\item{lambda}{lambda sequence.}

\item{algorithm.config}{the algorithm configuration to be used.}
}
\value{
\item{beta}{the fitted parameters -- the list \eqn{\hat\beta(\lambda(1)), \dots, \hat\beta(\lambda(d))} of length \code{length(return)}.
With each entry of list holding the fitted parameters, that is matrices of size \eqn{K\times p} (if \code{intercept = TRUE} matrices of size \eqn{K\times (p+1)})}
\item{loss}{the values of the loss function.}
\item{objective}{the values of the objective function (i.e. loss + penalty).}
\item{lambda}{the lambda values used.}
}
\description{
For a linear multiple output model with \eqn{p} features (covariates) dived into \eqn{m} groups using sparse group lasso.
}
\details{
This function computes a sequence of minimizers (one for each lambda given in the \code{lambda} argument) of
\deqn{\frac{1}{N}\|Y-X\beta\|_F^2 + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
where \eqn{\|\cdot\|_F} is the frobenius norm.
The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of features.
The group weights are denoted by \eqn{\gamma \in [0,\infty)^m} and the parameter weights by \eqn{\xi \in [0,\infty)^n}.
}
\author{
Martin Vincent
}

