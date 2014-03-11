#
#     Description of this R script:
#     R interface for linear sparse group lasso rutines.
#
#     Intended for use with R.
#     Copyright (C) 2013 Martin Vincent
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @title Fit a simple linear multi-response model using sparse group lasso 
#' 
#' @description
#' Fit a simple linear multi-response model with \eqn{p} covariates dived into groups using sparse group lasso.
#' 
#' @details
#' This function computes a sequence of minimizers (one for each lambda given in the \code{lambda} argument) of
#' \deqn{\|Y-X\beta\|_F^2 + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\|\cdot\|_F} is the frobenius norm.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of covariates.
#' The group weights are denoted by \eqn{\gamma \in [0,\infty)^m} and the parameter weights by \eqn{\xi = (\xi^{(1)},\dots, \xi^{(m)}) \in [0,\infty)^n}
#' with \eqn{\xi^{(1)}\in [0,\infty)^{n_1},\dots, \xi^{(m)} \in [0,\infty)^{n_m}}.
#' @param X design matrix, matrix of size \eqn{N \times p}.
#' @param Y response matrix, matrix of size \eqn{N \times K}.
#' @param weights, sample weights a vector of size \eqn{N}.
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda lambda sequence. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{beta}{the fitted parameters -- a list of length \code{length(return)} with each entry a matrix of size \eqn{q\times (p+1)} holding the fitted parameters.}
#' \item{loss}{the values of the loss function.}
#' \item{objective}{the values of the objective function (i.e. loss + penalty).}
#' \item{lambda}{the lambda values used.}
#' @author Martin Vincent
#' @useDynLib lsgl .registration=TRUE
#' @export
lsgl.multi <- function(X, Y, weights = rep(1/nrow(X), nrow(X)), covariateGrouping = factor(1:ncol(X)), groupWeights = NULL, parameterWeights =  NULL, alpha = 0.5, lambda, algorithm.config = sgl.standard.config) 
{
	
	#TODO dimension checks
	
	# Construct new weights
	weights <- rep(weights, ncol(Y))
	
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(ncol(Y)*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = ncol(Y), ncol = ncol(X))
	}
	
	# Construct design matrix
	X <- do.call(rbind, replicate(ncol(Y), X, simplify = FALSE))
	
	# Construct sample grouping
	G <- unlist(lapply(1:ncol(Y), function(g) rep(g, nrow(Y))))
	
	# Collapse Y to obtain response vector
	Y <- as.vector(Y)
	
	#Run lsgl
	res <- lsgl(x = X, y = Y, 
			weights = weights, 
			sampleGrouping = G,
			covariateGrouping = covariateGrouping, 
			groupWeights = groupWeights, 
			parameterWeights =  parameterWeights,
			alpha = alpha, 
			lambda = lambda,
			algorithm.config = algorithm.config)
	
	return(res)
	
}

#' Compute a lambda sequence for the regularization path
#' 
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param X design matrix, matrix of size \eqn{N \times p}.
#' @param Y response matrix, matrix of size \eqn{N \times K}.
#' @param weights, sample weights a vector of size \eqn{N}.
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence
#' @param lambda.min the smallest lambda value in the computed sequence. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return a vector of length \code{d} containing the compute lambda sequence.
#' @author Martin Vincent
#' @useDynLib lsgl .registration=TRUE
#' @export
lsgl.multi.lambda <- function(X, Y, weights = rep(1/(nrow(X)), nrow(X)), covariateGrouping = factor(1:ncol(X)), groupWeights = NULL, parameterWeights =  NULL, alpha = 0.5, d = 100L, lambda.min, algorithm.config = sgl.standard.config) 
{
	#TODO dimension checks

	# Construct new weights
	weights <- rep(weights, ncol(Y))
	
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(ncol(Y)*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = ncol(Y), ncol = ncol(X))
	}
	
	# Construct design matrix
	X <- do.call(rbind, replicate(ncol(Y), X, simplify = FALSE))
	
	# Construct sample grouping
	G <- unlist(lapply(1:ncol(Y), function(g) rep(g, nrow(Y))))
	
	# Collapse Y to obtain response vector
	Y <- as.vector(Y)
	
	#Run lsgl
	res <- lsgl.lambda(x = X, y = Y, 
			weights = weights, 
			sampleGrouping = G,
			covariateGrouping = covariateGrouping, 
			groupWeights = groupWeights, 
			parameterWeights =  parameterWeights,
			alpha = alpha, 
			d = d, 
			lambda.min = lambda.min,
			algorithm.config = algorithm.config)
	
	return(res)	
}


#' Cross validation using multiple possessors 
#' 
#' @param X design matrix, matrix of size \eqn{N \times p}.
#' @param Y response matrix, matrix of size \eqn{N \times K}.
#' @param weights, sample weights a vector of size \eqn{N}.
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda lambda sequence. 
#' @param fold
#' @param cv.indices
#' @param max.threads
#' @param algorithm.config the algorithm configuration to be used. 
#' @return sgl object.
#' @author Martin Vincent
#' @useDynLib lsgl .registration=TRUE
#' @export
lsgl.multi.cv <- function(X, Y, weights = rep(1/nrow(X), nrow(X)), 
		covariateGrouping = factor(1:ncol(X)), 
		groupWeights = NULL,
		parameterWeights =  NULL, 
		alpha = 0.5, lambda, fold = 10L, cv.indices = list(), max.threads = 2L,
		algorithm.config = sgl.standard.config) 
{
	#TODO dimension checks
		
	# Construct new weights
	weights <- rep(weights, ncol(Y))
	
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(ncol(Y)*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = ncol(Y), ncol = ncol(X))
	}
	
	# Construct design matrix
	X <- do.call(rbind, replicate(ncol(Y), X, simplify = FALSE))
	
	# Construct sample grouping
	G <- unlist(lapply(1:ncol(Y), function(g) rep(g, nrow(Y))))
	
	# Collapse Y to obtain response vector
	Y <- as.vector(Y)
	
	#Run lsgl
	res <- lsgl.cv(x = X, y = Y, 
			weights = weights, 
			sampleGrouping = G,
			covariateGrouping = covariateGrouping, 
			groupWeights = groupWeights, 
			parameterWeights =  parameterWeights,
			alpha = alpha, 
			lambda = lambda,
			fold = fold,
			cv.indices = cv.indices,
			max.threads = max.threads,
			algorithm.config = algorithm.config)
	
	return(res)
	
}
