# TODO: Add comment
# 
# Author: martin
###############################################################################

rearrange <- function(data, covariate.order, ...) UseMethod("rearrange")


#' Rearrange sgl.data
#' 
#' @param data 
#' @param covarite.order 
#' @author Martin Vincent
#' @export
rearrange.sgldata <- function(data, covariate.order, ...) 
{
	
	data$x <- data$x[,covariate.order]
	data$covariate.names <- data$covarite.names[covariate.order]
	
	return(data)
}

#' Create sgl data
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param y responses, vector of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param sampleGrouping grouping of samples, a factor of length \eqn{N}. Default is no grouping (NULL), that is all samples is the same group
#' @author Martin Vincent
#' @export
create.sgldata <- function(x, y, weights = rep(1/nrow(x), nrow(x)), sampleGrouping = NULL) {
	
	#TODO dim checks
	
	data <- list()
	data$X <- as.matrix(x)
	data$Y <- as.numeric(y)
	data$W <- as.numeric(weights)
	
	# sample grouping
	
	if(is.null(sampleGrouping)) {
		sampleGrouping <- rep(1, nrow(x))
	}
	
	sampleGrouping <- factor(sampleGrouping)
	data$G <- as.integer(factor(sampleGrouping))-1L
	
	# dimensions
	data$n.covariate <- ncol(x)
	data$n.groups <- length(levels(sampleGrouping))
	
	# names
	data$sample.names <- rownames(x)
	data$covariate.names <- colnames(x)
	data$group.names <- levels(sampleGrouping)
	
	# Is X sparse
	data$sparseX <- is(data$X, "sparseMatrix")
	
	class(data) <- "sgldata"
	return(data)
}

#' Prepare sgl function arguments 
#' 
#' @param data 
#' @param covariateGrouping 
#' @param groupWeights 
#' @param parameterWeights 
#' @param alpha 
#' @author Martin Vincent
#' @export
prepare.args <- function(data, covariateGrouping, groupWeights, parameterWeights, alpha) {
	
	#If Lasso then ignore grouping
	if(alpha == 1) {
		covariateGrouping <- factor(1:data$n.covariate)
		groupWeights <- rep(1, data$n.covariate)
	}
	
	ncov <- data$n.covarites
	ngrp <- data$n.groups
		
	group.order <- order(covariateGrouping)
	
	#Reorder data
	data <- rearrange(data, group.order)
	
	parameterWeights <- parameterWeights[,group.order, drop = FALSE]
	
	#Compute block dim
	block.dim <- ngrp*as.integer(table(covariateGrouping))
	
	#args list
	args <- list()
	args$block.dim <- block.dim
	args$groupWeights <- groupWeights
	args$parameterWeights <- parameterWeights
	args$alpha <- alpha
	args$data <- data
	args$group.order <- group.order
	
	return(args)
}