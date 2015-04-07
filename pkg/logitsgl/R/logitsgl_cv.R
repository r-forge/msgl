#
#     Description of this R script:
#     R interface for linear multiple output sparse group lasso routines.
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
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

#' @title Linear multiple output cross validation using multiple possessors 
#' 
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y response matrix, matrix of size \eqn{N \times K}.
#' @param intercept should the model include intercept parameters
#' @param grouping grouping of features, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the feature. 
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than 1 and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' @param cv.indices a list of indices of a cross validation splitting.
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads the maximal number of threads to be used.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return
#' \item{Yhat}{the cross validation estimated response matrix}
#' \item{Y.true}{the true response matrix, this is equal to the argument \code{y}}
#' \item{cv.indices}{the cross validation splitting used}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models.}
#' @author Martin Vincent
#' @useDynLib logitsgl .registration=TRUE
#' @export
logitsgl.cv <- function(x, y, intercept = TRUE,
		grouping = factor(1:ncol(x)), 
		groupWeights = c(sqrt(ncol(y)*table(grouping))),
		parameterWeights =  matrix(1, nrow = ncol(y), ncol = ncol(x)), 
		alpha = 1, lambda, fold = 10L, cv.indices = list(), max.threads = 2L,
		algorithm.config = logitsgl.standard.config) 
{
	# Get call
	cl <- match.call()
	
	# cast
	grouping <- factor(grouping)
	
	# add intercept
	if(intercept) {
		x <- cBind(Intercept = rep(1, nrow(x)), x)
		groupWeights <- c(0, groupWeights)
		parameterWeights <- cbind(rep(0, ncol(y)), parameterWeights)
		grouping <- factor(c("Intercept", as.character(grouping)), levels = c("Intercept", levels(grouping)))
	}
	
	# create data
	group.names <- if(is.null(colnames(y))) 1:ncol(y) else colnames(y)
	data <- create.sgldata(x, y, group.names = group.names)
	
	# call SglOptimizer function
	callsym <- paste("logitsgl_", if(data$sparseX) "xs_" else "xd_", if(data$sparseY) "ys" else "yd", sep = "")
	res <- sgl_cv(callsym, "logitsgl", data, grouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, algorithm.config)
	
	# Add true response
	res$Y.true <- y
	
	# Responses
	res$P <- lapply(res$responses$prob, t)
	res$link <- lapply(res$responses$link, t)
	res$Yhat <- lapply(res$responses$classes, t)
	res$responses <- NULL
	
	#TODO response dimnames
		
	res$logitsgl_version <- packageVersion("logitsgl")
	res$intercept <- intercept
	res$call <- cl
	
	class(res) <- "logitsgl"
	return(res)
}
