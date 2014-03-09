#
#     Description of this R script:
#     R interface to sgl-predict
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

#' Sgl predict
#' 
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param object a sgl object containing a list of estimated models.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param ... not used.
#' @return sgl object content will depend on the C++ response class.
#' @author Martin Vincent
#' @export
sgl_predict <- function(module_name, PACKAGE, object, data, ...) {

	if("beta" %in% names(object)) {

		beta <- lapply(X = object$beta, FUN = function(m) as(m, "CsparseMatrix"))
		beta <- lapply(X = beta, FUN = function(m) list(dim(m), m@p, m@i, m@x))
		
		call_sym <- paste(module_name, "sgl_predict", sep="_")
		res <- .Call(call_sym, PACKAGE = PACKAGE, data, beta)
		
	} else  {
                stop("No models found -- missing beta")
	}

	if(!is.null(data$sample.names)) {
		# Set sample names
		res$responses <- lapply(res$responses, function(x) .set_sample_names(x, data$sample.names))
	}
	
	res$sglOptim_version <- packageVersion("sglOptim")
	class(res) <- "sgl"
	
	return(res)
}
