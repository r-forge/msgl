#
#     Description of this R script:
#     R interface for linear multi-response sparse group lasso routines.
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

#' @title Predict
#'
#' @description 
#' Compute the predicted response matrix for a new data set.
#'
#' @param object an object of class logitsgl, produced with \code{logitsgl}.
#' @param x a data matrix of size \eqn{N_\textrm{new} \times p}.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param ... ignored.
#' @return
#' \item{Yhat}{the predicted response matrix (of size \eqn{N_\textrm{new} \times K})}
#' @author Martin Vincent
#' @method predict logitsgl
#' @export
#' @useDynLib logitsgl .registration=TRUE
predict.logitsgl <- function(object, x, sparse.data = is(x, "sparseMatrix"), ...) 
{
	# Get call
	cl <- match.call()
	
	if(object$intercept){
		# add intercept
		x <- cBind(Intercept = rep(1, nrow(x)), x)
	}	
	
	data <- list()
	
	if(sparse.data) {
		
		x <- as(x, "CsparseMatrix")
		data$X <- list(dim(x), x@p, x@i, x@x)
		
		res <- sgl_predict("logitsgl_xs_yd", "logitsgl", object, data)
		
	} else {
		
		data$X <- as.matrix(x)
		
		res <- sgl_predict("logitsgl_xd_yd", "logitsgl", object, data)
		
	}
	
	#Responses
	res$P <- lapply(res$responses$prob[,1], t)
	res$link <- lapply(res$responses$link, t)
	res$Yhat <- lapply(res$responses$classes, t)
	res$responses <- NULL
	
	#TODO response dimnames
	
	res$call <- cl
	
	return(res)
}
