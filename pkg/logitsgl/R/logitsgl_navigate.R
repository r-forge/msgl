#
#     Description of this R script:
#     R interface for linear multi-response models using sparse group lasso.
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

#TODO loss function
#' @title Compute error rates 
#' 
#' @description
#' Compute and return the error
#' TODO
#' for each model.
#' 
#' @param object a logitsgl object
#' @param data a design matrix (the \eqn{X} matrix)
#' @param response redirected to \code{y}
#' @param y a matrix of the true responses (the \eqn{Y} matrix)
#' @param ... ignored
#' @return a vector of error rates
#' 
#' @author Martin Vincent
#' @method Err logitsgl
#' @import sglOptim
#' @export
Err.logitsgl <- function(object, data = NULL, response = object$Y.true, Y = response, type = "rate", ... ) {

	if(type=="rate") {
		return(compute_error(object, data = data, response.name = "Yhat", response = Y, loss = function(x,y) mean(x != y)))
	}
	
	if(type=="count") {
		return(compute_error(object, data = data, response.name = "Yhat", response = Y, loss = function(x,y) sum(x != y)))
	}
	
	if(type=="loglike") {
		loss <- function(x,y) -mean(y*log(x)+(1-y)*log(1-x))
		return(compute_error(object, data = data, response.name = "P", response = Y, loss = loss))
	}
	
	stop("Unknown type")
	
}


#' @title Nonzero features
#' 
#' @description
#' Extracts the nonzero features for each model.
#'
#' @param object a logitsgl object
#' @param ... ignored
#' @return a list of of length \code{nmod(x)} containing the nonzero features (that is nonzero columns of the beta matrices)
#' @author Martin Vincent
#' @method features logitsgl
#' @import sglOptim
#' @export
features.logitsgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(features(object))
}

#' @title Nonzero parameters
#' 
#' @description
#' Extracts the nonzero parameters for each model.
#'
#' @param object a logitsgl object
#' @param ... ignored
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#' @author Martin Vincent
#' @method parameters logitsgl
#' @import sglOptim
#' @export
parameters.logitsgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(parameters(object))
}

#' @title Returns the number of models in a lsgl object
#'
#' @param object a logitsgl object
#' @param ... ignored
#' @return the number of models in \code{object}
#' @author Martin Vincent
#' @method nmod logitsgl
#' @import sglOptim
#' @export
nmod.logitsgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(nmod(object, ...))
}

#' @title Exstract the fitted models 
#' 
#' @description
#' Returns the fitted models, that is the estimated \eqn{\beta} matrices.
#' 
#' @param object a logitsgl object 
#' @param index indices of the models to be returned
#' @param ... ignored
#' @return a list of \eqn{\beta} matrices.
#' 
#' @author Martin Vincent
#' @method models logitsgl
#' @import sglOptim
#' @export
models.logitsgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	return(models(object, ...))
}

#' @title Extract nonzero coefficients 
#'
#' @param object a logitsgl object
#' @param index indices of the models
#' @param ... ignored
#' @return a list of length \code{length(index)} with nonzero coefficients of the models
#'
#' @author Martin Vincent
#' @method coef logitsgl
#' @import sglOptim
#' @export
coef.logitsgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	return(coef(object, index = index, ...))
}


#' Print function for logitsgl
#'
#' This function will print some general information about the lsgl object
#'  
#' @param x logitsgl object
#' @param ... ignored
#' 
#' @method print logitsgl
#' @author Martin Vincent
#' @import sglOptim
#' @export
print.logitsgl <- function(x, ...) {
	sgl_print(x)
}
