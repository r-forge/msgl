# TODO: Add comment
# 
# Author: martin
###############################################################################

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

	res <- list()
	
	if("beta" %in% names(object)) {

		beta <- lapply(X = object$beta, FUN = function(m) as(m, "CsparseMatrix"))
		beta <- lapply(X = beta, FUN = function(m) list(dim(m), m@p, m@i, m@x))
		
		call_sym <- paste(module_name, "sgl_predict", sep="_")
		res <- .Call(call_sym, PACKAGE = PACKAGE, data, beta)
		
	} else  {
		stop("No models found - missing beta")
	}
	
	class(res) <- "sgl"
	return(res)
}