% Generated by roxygen2 (4.0.1): do not edit by hand
\name{msgl.lambda.seq}
\alias{msgl.lambda.seq}
\title{Computes a lambda sequence for the regularization path}
\usage{
msgl.lambda.seq(x, classes, sampleWeights = rep(1/length(classes),
  length(classes)), grouping = NULL, groupWeights = NULL,
  parameterWeights = NULL, alpha = 0.5, d = 100L, standardize = TRUE,
  lambda.min, intercept = TRUE, sparse.data = is(x, "sparseMatrix"),
  algorithm.config = sgl.standard.config)
}
\arguments{
\item{x}{design matrix, matrix of size \eqn{N \times p}.}

\item{classes}{classes, factor of length \eqn{N}.}

\item{sampleWeights}{sample weights, a vector of length \eqn{N}.}

\item{grouping}{grouping of features, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate.}

\item{groupWeights}{the group weights, a vector of length \eqn{m+1} (the number of groups).
The first element of the vector is the intercept weight.
If \code{groupWeights = NULL} default weights will be used.
Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of features in the group}}} for all other weights.}

\item{parameterWeights}{a matrix of size \eqn{K \times (p+1)}.
The first column of the matrix is the intercept weights.
Default weights are is 0 for the intercept weights and 1 for all other weights.}

\item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.}

\item{d}{the length of lambda sequence}

\item{standardize}{if TRUE the features are standardize before fitting the model. The model parameters are returned in the original scale.}

\item{lambda.min}{the smallest lambda value in the computed sequence.}

\item{intercept}{should the model include intercept parameters}

\item{sparse.data}{if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.}

\item{algorithm.config}{the algorithm configuration to be used.}
}
\value{
a vector of length \code{d} containing the computed lambda sequence.
}
\description{
Computes a decreasing lambda sequence of length \code{d}.
The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100, lambda.min = 0.01)
}
\author{
Martin Vincent
}

