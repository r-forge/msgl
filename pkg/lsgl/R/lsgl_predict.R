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

#' Predict
#' 
#' @param object 
#' @param x 
#' @param ... 
#' @return
#' \item{Yhat}{the linear predictors -- a list of length \code{length(fit$beta)} one item for each model, with each item a matrix of size \eqn{N_\textrm{new} \times K} containing the linear predictors.}
#' @author Martin Vincent
#' @export
predict.lsgl <- function(object, x, ...) 
{
	
	if(object$intercept){
		# add intercept
		x <- cBind(Intercept = rep(1, nrow(x)), x)
	}	
	
	object$beta <- lapply(object$beta, t)
	
	data <- list()
	
	if(is(x, "sparseMatrix")) {
		
		x <- as(x, "CsparseMatrix")
		data$X <- list(dim(x), x@p, x@i, x@x)
		
		res <- sgl_predict("lsgl_sparse", "lsgl", object, data)
		
	} else {
		
		data$X <- as.matrix(x)
		
		res <- sgl_predict("lsgl_dense", "lsgl", object, data)
		
	}
	
	#Responses
	
	res$Yhat <- lapply(res$responses$link, t)
	res$responses <- NULL
	
	return(res)
}
