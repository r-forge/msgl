# TODO: Add comment
# 
# Author: martin
###############################################################################
#' Sparse group lasso cross validation using multiple possessors 
#' 
#' @param call_sym reference to objective specific C++ routines
#' @param data
#' @param covariateGrouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than \eqn{1} and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' If \code{fold}\eqn{\le}\code{max(table(classes))} then the data will be split into \code{fold} disjoint subsets keeping the ration of classes approximately equal.
#' Otherwise the data will be split into \code{fold} disjoint subsets without keeping the ration fixed.
#' @param cv.indices a list of indices of a cross validation splitting. 
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads the maximal number of threads to be used
#' @param seed the seed used for generating the random cross validation splitting, only used if \code{fold}\eqn{\le}\code{max(table(classes))}. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @author Martin Vincent
#' @export
sgl_cv <- function(module_name, PACKAGE, data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold = 2L, cv.indices = list(), max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
	args <- prepare.args(data, covariateGrouping, groupWeights, parameterWeights, alpha)
	
	
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
