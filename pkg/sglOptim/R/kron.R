
#' Create Kronecker product representation 
#' 
#' @param ... matrices in the product
#' @return Kronecker product representation 
#' 
#' @author Martin Vincent
#' @export
kron <- function(...) {
	
	x <- lapply(as.list(match.call())[-1], function(q) eval(q))
	
	if(length(x) == 0) {
		stop("argument must consist of matrices")
	}
	
	if(!all(sapply(x, function(M) is(M, "matrix")))) {
		stop("argument must consist of matrices (sparse matrices not supported)")
	}
	
	if(!(length(x) %in% c(2,3))) {
		stop("only dual and triple Kronecker products are supported")
	}
	
	class(x) <- "kron"

	return(x)
}

#' Convert a list of matrices to a kron object
#' @param x a list of matrices
#' @return a kron object
#' 
#' @author martin
#' @export
as.kron <- function(x) {
	
	if(!is.list(x)) {
		stop("argument must be list of matrices")
	}
	
	if(!all(sapply(x, function(M) is(M, "matrix")))) {
		stop("argument must consist of matrices (sparse matrices not supported)")
	}
	
	class(x) <- "kron"

	return(x)
}

#' dimension of kron object
#' @param x kron object
#' @return the dimension
#' 
#' @author Martin Vincent
#' @method dim kron
#' @export
dim.kron <- function(x) {
	return(c(prod(sapply(x, nrow)),prod(sapply(x, ncol))))
}
