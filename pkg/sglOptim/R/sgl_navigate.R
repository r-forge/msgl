#
#     Description of this R script:
#     Routines for navigating sgl objects
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

# S3 functions:

#' Generic function for extracting nonzero features (or groups) 
#' @param x an object
#' @param ... additional paramters (optional)
#' @return a list of nonzero features
#' 
#' @author Martin Vincent
#' @export
features <- function(x, ...) UseMethod("features")

#' Generic function for extracting nonzero parameters
#' @param x an object
#' @param ... additional paramters (optional)
#' @return a list of nonzero paramters
#' 
#' @author Martin Vincent
#' @export
parameters <- function(x, ...) UseMethod("parameters")


#' Generic function for couting the number of models
#' @param x an object
#' @param ... additional paramters (optional)
#' @return the number of models contained in \code{x}
#' 
#' @author Martin Vincent
#' @export
nmod <- function(x, ...) UseMethod("nmod")

#' Generic function for extracting the fitted models 
#' @param x an object
#' @param index a vector of indices of the models to be returned
#' @param ... additional paramters (optional)
#' @return a list of length \code{length(index)} containing the models
#' 
#' @author Martin Vincent
#' @export
models <- function(x, index, ...) UseMethod("models")


#' Generic function for extracting the nonzero coefficients 
#' @param x an object
#' @param index a vector of indices of the models 
#' @param ... additional paramters (optional)
#' @return a list of length \code{length(index)} containing the nonzero paramters of the models
#' 
#' @author Martin Vincent
#' @export
coef <- function(x, index, ...) UseMethod("coef")


#' todo
#' @param x 
#' @param ... 
#' @return a list of nonzero features (that is nonzero colums of the beta matrices)
#' 
#' @author martin
#' @method features sgl
#' @S3method features sgl
#' @export
features.sgl <- function(x, ...) {
	
	if(is.null(x$beta)) {
		stop("object contains no models")
	}
	
	if(is.null(colnames(x$beta[[1]]))) {
		res <- lapply(x$beta, function(beta) which(colSums(beta != 0) != 0))
	} else {
		res <- lapply(x$beta, function(beta) colnames(beta)[colSums(beta != 0) != 0])
	}
	
	return(res)
}

#' todo
#' @param x 
#' @param ... 
#' @return todo
#' 
#' @author martin
#' @method parameters sgl
#' @S3method parameters sgl
#' @export
parameters.sgl <- function(x, ...) {
	
	if(is.null(x$beta)) {
		stop("object contains no models")
	}
	
	tmp <- features(x)
	res <- sapply(1:length(x$beta), function(i) x$beta[[i]][,tmp[[i]]] != 0)
		
	return(res)
}

#' todo
#' @param x 
#' @param ... 
#' @return todo
#' 
#' @author martin
#' @method nmod sgl
#' @S3method nmod sgl
#' @export
nmod.sgl <- function(x, ...) {
	
	if(is.null(x$beta)) {
		return (0)
	}
	
	return(length(x$beta))
}

#' Exstract the estimated models
#' 
#' @param x a msgl object 
#' @param index the models to be returned
#' @return a list of sparse matrices
#' 
#' @author Martin Vincent
#' @method models sgl
#' @S3method models sgl
#' @export
models.sgl <- function(x, index = 1:nmod(x), ...) {
	
	if(is.null(x$beta)) {
		stop("object contains no models")
	}
	
	return(x$beta[index])
}

#' todo
#' @param x 
#' @param index 
#' @param ... 
#' @return todo
#' 
#' @author martin
#' @method coef sgl
#' @S3method coef sgl
#' @export
coef.sgl <- function(x, index = 1:nmod(x), ...) {
	
	if(is.null(x$beta)) {
		stop("object contains no models")
	}
	
	return(lapply(x$beta[index], function(beta) beta[,colSums(beta != 0) != 0]))
}

