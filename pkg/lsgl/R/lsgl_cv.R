#
#     Description of this R script:
#     R interface for multinomial sparse group lasso rutines.
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

#' Cross validation using multiple possessors 
#' 
#' @param x
#' @param y
#' @param intercept
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambd lambda sequence. 
#' @param fold
#' @param cv.indices
#' @param max.threads
#' @param algorithm.config the algorithm configuration to be used. 
#' @return sgl object.
#' @author Martin Vincent
#' @useDynLib lsgl .registration=TRUE
#' @export
lsgl.cv <- function(x, y, intercept = TRUE,
		covariateGrouping = factor(1:ncol(x)), 
		groupWeights = c(sqrt(ncol(y)*table(covariateGrouping))),
		parameterWeights =  matrix(1, nrow = ncol(y), ncol = ncol(x)), 
		alpha = 0.5, lambda, fold = 10L, cv.indices = list(), max.threads = 2L,
		algorithm.config = sgl.standard.config) 
{
	# cast
	covariateGrouping <- factor(covariateGrouping)
	sampleGrouping <- factor(sampleGrouping)
	
	# add intercept
	if(intercept) {
		x <- cBind(Intercept = rep(1, nrow(x)), x)
		groupWeights <- c(0, groupWeights)
		parameterWeights <- cbind(rep(0, ncol(y)), parameterWeights)
		covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	}
	
	# create data
	data <- create.sgldata(x, y, sampleGrouping = 1:ncol(y))
	
	# call SglOptimizer function
	if(data$sparseX) {
		
		res <- sgl_cv("lsgl_sparse", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, algorithm.config)
		
	} else {
		
		res <- sgl_cv("lsgl_dense", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, algorithm.config)
		
	}
	
	res$intercept <- intercept
	
	class(res) <- "lsgl"
	return(res)
}
