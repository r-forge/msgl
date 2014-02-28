#
#     Description of this R script:
#     Routines for handling sgl-data
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


rearrange <- function(data, covariate.order, ...) UseMethod("rearrange")


#' Rearrange sgldata
#' 
#' @param data  sgldata object
#' @param covariate.order the new order of the covarites
#' @param ... not used
#' @return a sgl data object with the covariates reordered
#' @author Martin Vincent
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
#' @param weights sample weights, a vector of length \eqn{N}.
#' @param sampleGrouping grouping of samples, a factor of length \eqn{N}. Default is no grouping (NULL), that is all samples is the same group.
#' @param sparseX if TRUE \code{x} will be treated as sparse, if FALSE \code{x} will be treated as dens.
#' @author Martin Vincent
#' @export
create.sgldata <- function(x, y, weights = rep(1/nrow(x), nrow(x)), sampleGrouping = NULL, sparseX = is(x, "sparseMatrix")) {
	
	#TODO dim checks
	
	data <- list()
	
	# Is X sparse
	data$sparseX <- sparseX

	if(data$sparseX) {
		data$X <- as(x, "CsparseMatrix")
	} else {
		data$X <- as.matrix(x)
	}

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

	# sparse X format
	if(data$sparseX) {
		data$X <- list(dim(data$X), data$X@p, data$X@i, data$X@x)
	}
		
	class(data) <- "sgldata"
	return(data)
}

#' Prepare sgl function arguments 
#' 
#' @param data sgldata object
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @author Martin Vincent
prepare.args <- function(data, parameterGrouping, groupWeights, parameterWeights, alpha) {
	
	#If Lasso then ignore grouping
	if(alpha == 1) {
		parameterGrouping <- factor(1:data$n.covariate)
		groupWeights <- rep(1, data$n.covariate)
	}
	
	ncov <- data$n.covarites
	ngrp <- data$n.groups
		
	group.order <- order(parameterGrouping)
	
	#Reorder data
	data <- rearrange(data, group.order)
	
	parameterWeights <- parameterWeights[,group.order, drop = FALSE]
	
	#Compute block dim
	block.dim <- ngrp*as.integer(table(parameterGrouping))
	
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
