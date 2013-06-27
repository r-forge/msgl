# TODO: Add comment
# 
# Author: martin
###############################################################################

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

#' Fit
#' 
#' @param x
#' @param y
#' @param weights
#' @param sampleGrouping
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambd lambda sequence. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return a vector of length \code{d} containing the compute lambda sequence.
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

#' Cross validation using multiple possessors 
#' 
#' @param x
#' @param y
#' @param weights
#' @param sampleGrouping
#' @param covariateGrouping grouping of covariates, a factor or vector of length \eqn{p}. Each element of the factor/vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambd lambda sequence. 
#' @param fold
#' @param cv.indices
#' @param max.threads
#' @param seed
#' @param algorithm.config the algorithm configuration to be used. 
#' @return sgl object.
#' @author Martin Vincent
#' @useDynLib lsgl .registration=TRUE
#' @export
lsgl.cv <- function(x, y, weights = rep(1/nrow(x), nrow(x)), 
		sampleGrouping = factor(rep(1, nrow(x))), 
		covariateGrouping = factor(1:ncol(x)), 
		groupWeights = c(sqrt(length(levels(sampleGrouping))*table(covariateGrouping))),
		parameterWeights =  matrix(1, nrow = length(levels(sampleGrouping)), ncol = ncol(x)), 
		alpha = 0.5, lambda, fold = 10L, cv.indices = list(), max.threads = 2L, seed = 331L, 
		algorithm.config = sgl.standard.config) 
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
		
		res <- sgl_cv("lsgl_sparse", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, seed, algorithm.config)
		
	} else {
		
		res <- sgl_cv("lsgl_dense", "lsgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, seed, algorithm.config)
	
	}
	
	class(res) <- "lsgl"
	return(res)
}

#' Predict
#' 
#' @param object 
#' @param x 
#' @param ... 
#' @return sgl object
#' @author Martin Vincent
#' @export
predict.lsgl <- function(object, x, ...) 
{
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	
	data <- list()
	
	if(is(x, "sparseMatrix")) {
		
		x <- as(x, "CsparseMatrix")
		data$X <- list(dim(x), x@p, x@i, x@x)

		res <- sgl_predict("lsgl_sparse", "lsgl", object, data)
		
	} else {
		
		data$X <- as.matrix(x)

		res <- sgl_predict("lsgl_dense", "lsgl", object, data)
		
	}
	
	
	return(res)
}