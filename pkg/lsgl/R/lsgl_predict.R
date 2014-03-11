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
