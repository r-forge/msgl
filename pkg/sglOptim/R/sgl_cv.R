#
#     Description of this R script:
#     TODO
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

#' Generic sparse group lasso cross validation using multiple possessors 
#' 
#' 
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than \eqn{1} and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' If \code{fold}\eqn{\le}\code{max(table(classes))} then the data will be split into \code{fold} disjoint subsets keeping the ration of classes approximately equal.
#' Otherwise the data will be split into \code{fold} disjoint subsets without keeping the ration fixed.
#' @param cv.indices a list of indices of a cross validation splitting. 
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads the maximal number of threads to be used.
#' @param seed the seed used for generating the random cross validation splitting, only used if \code{fold}\eqn{\le}\code{max(table(classes))}. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return sgl object content will depend on the C++ response class.
#' @author Martin Vincent
#' @export
sgl_cv <- function(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, lambda, fold = 2L, cv.indices = list(), max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
	args <- prepare.args(data, parameterGrouping, groupWeights, parameterWeights, alpha)
	
	
	if(length(cv.indices) == 0) {
		
		use.cv.indices <- FALSE
		
		# Check fold
		if(fold < 2) {
			stop("fold must be equal to or larger than 2")
		}
		
		if(fold > length(data$G)) {
			stop("fold must be equal to or less than the number of samples")
		}
		
		if(fold > max(table(data$G))) {
			# use random sample indices
			use.cv.indices <- TRUE
			cv.indices <- split(sample(0:(length(data$G))-1L), 1:fold)
		}
		
	} else {
		
		cv.indices <- lapply(cv.indices, function(x) as.integer(x-1))
		use.cv.indices <- TRUE
	}
	
	call_sym <- paste(module_name, "sgl_cv", sep="_")
	res <- .Call(call_sym, PACKAGE = PACKAGE, args$data, args$block.dim, args$groupWeights, args$parameterWeights, args$alpha, lambda, fold, cv.indices, use.cv.indices, max.threads, seed, algorithm.config)				
		
	class(res) <- "sgl"
	return(res)
}
