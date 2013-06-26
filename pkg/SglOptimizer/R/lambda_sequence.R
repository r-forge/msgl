# TODO: Add comment
# 
# Author: martin
###############################################################################



#' Generic routine for computing a lambda sequence for the regularization path
#' 
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param call_sym reference to objective specific C++ routines
#' @param PACKAGE
#' @param data
#' @param sampleGrouping
#' @param covariateGrouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence
#' @param lambda.min the smallest lambda value in the computed sequence. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return a vector of length \code{d} containing the compute lambda sequence.
#' @author Martin Vincent
#' @export
#' @import Matrix
sgl_lambda_sequence <- function(module_name, PACKAGE, data, covariateGrouping, groupWeights, parameterWeights, alpha = 0.5, d = 100L, lambda.min, algorithm.config = sgl.standard.config) {
	
	args <- prepare.args(data, covariateGrouping, groupWeights, parameterWeights, alpha)
	
	call_sym <- paste(module_name, "sgl_lambda", sep="_")
		res <- .Call(call_sym, PACKAGE = PACKAGE, args$data, args$block.dim, args$groupWeights, args$parameterWeights, args$alpha, d, lambda.min, algorithm.config)
	
	return(res)
}
