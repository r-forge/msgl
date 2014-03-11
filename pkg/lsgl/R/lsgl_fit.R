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

#' @title Fit a linear multi-response model using sparse group lasso 
#' 
#' @description
#' For a linear multi-response model with \eqn{p} covariates dived into \eqn{m} groups using sparse group lasso.
#' 
#' @details
#' This function computes a sequence of minimizers (one for each lambda given in the \code{lambda} argument) of
#' \deqn{R(\beta) + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{R} is given below.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of covariates.
#' The group weights are denoted by \eqn{\gamma \in [0,\infty)^m} and the parameter weights by \eqn{\xi = (\xi^{(1)},\dots, \xi^{(m)}) \in [0,\infty)^n}
#' with \eqn{\xi^{(1)}\in [0,\infty)^{n_1},\dots, \xi^{(m)} \in [0,\infty)^{n_m}}.
#'
#' \deqn{R(\beta) = \sum_{i=1}^N w_i \left(y_i - \beta_{0g_i}\sum_{j=1}^p \beta_{g_ij}x_j\right)}
#' where \eqn{w} is the sample weights and \eqn{g} the sample grouping.
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response, vector of size \eqn{N}.
#' @param weights, sample weights a vector of size \eqn{N}.
#' @param sampleGrouping, vector of size \eqn{N}
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{beta}{the fitted parameters -- a list of length \code{length(lambda)} with each entry a matrix of size \eqn{K\times (p+1)}.}
#' \item{loss}{the values of the loss function.}
#' \item{objective}{the values of the objective function (i.e. loss + penalty).}
#' \item{lambda}{the lambda values used.}
#' @author Martin Vincent
#' @useDynLib lsgl .registration=TRUE
#' @export
lsgl <- function(x, y, weights = rep(1/nrow(x), nrow(x)), sampleGrouping = factor(rep(1, nrow(x))), covariateGrouping = factor(1:ncol(x)), groupWeights = c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping))), parameterWeights =  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x)), alpha = 0.5, lambda, algorithm.config = sgl.standard.config) 
{
	# cast
	covariateGrouping <- factor(covariateGrouping)
	sampleGrouping <- factor(sampleGrouping)
	
	# add intercept
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(sampleGrouping))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	
	# create data
	data <- create.sgldata(x, y, weights, sampleGrouping)
	
	# call SglOptimizer function
	if(data$sparseX) {
		res <- sgl_fit("lsgl_sparse", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)
	} else {
		res <- sgl_fit("lsgl_dense", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)
	}
	
	class(res) <- "lsgl"
	return(res)
}

#' Compute a lambda sequence for the regularization path
#' 
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param x
#' @param y
#' @param weights
#' @param sampleGrouping
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
lsgl.lambda <- function(x, y, weights = rep(1/nrow(x), nrow(x)), sampleGrouping = factor(rep(1, nrow(x))), covariateGrouping = factor(1:ncol(x)), groupWeights = c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping))), parameterWeights =  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x)), alpha = 0.5, d = 100L, lambda.min, algorithm.config = sgl.standard.config) 
{
	# cast
	covariateGrouping <- factor(covariateGrouping)
	sampleGrouping <- factor(sampleGrouping)
	
	# add intercept
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(sampleGrouping))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	
	# create data
	data <- create.sgldata(x, y, weights, sampleGrouping)
	
	# call SglOptimizer function
	if(data$sparseX) {
		lambda <- sgl_lambda_sequence("lsgl_sparse", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
	} else {
		lambda <- sgl_lambda_sequence("lsgl_dense", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
	}
	
	return(lambda)
}

