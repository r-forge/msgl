# TODO: Add comment
# 
# Author: martin
###############################################################################

#' Fit a sparse group lasso regularization path. 
#'
#' @param call_sym reference to objective specific C++ routines
#' @param data
#' @param covariateGrouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param return the indices of lambda values for which to return a the fitted parameters.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{beta}{the fitted parameters -- a list of length \code{length(lambda)} with each entry a matrix of size \eqn{K\times (p+1)} holding the fitted parameters}
#' \item{loss}{the values of the loss function}
#' \item{objective}{the values of the objective function (i.e. loss + penalty)}
#' \item{lambda}{the lambda values used}
#' @author Martin Vincent
#' @export
#' @import Matrix
sgl_fit <- function(module_name, PACKAGE, data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config = sgl.standard.config) {
	
	args <- prepare.args(data, covariateGrouping, groupWeights, parameterWeights, alpha)
	
	#TODO check that return is valid
	return <- as.integer(sort(unique(return))) - 1L
	
	call_sym <- paste(module_name, "sgl_fit", sep="_")
	res <- .Call(call_sym, PACKAGE = PACKAGE, args$data, args$block.dim, args$groupWeights, args$parameterWeights, args$alpha, lambda, return, algorithm.config)
	
	# Dim names
	covariate.names <- args$data$covariate.names
	group.names <- args$data$group.names
	
	# Create R sparse matrix
	res$beta <- lapply(1:length(res$beta), function(i) sparseMatrix(p = res$beta[[i]][[2]], i = res$beta[[i]][[3]], x = res$beta[[i]][[4]], dims = res$beta[[i]][[1]], dimnames = list(group.names, covariate.names), index1 = FALSE))
	
	# Restore org order
	res$beta <- lapply(res$beta, function(beta.matrix) beta.matrix[, order(args$group.order)])
		
	class(res) <- "sgl"
	return(res)
}