#
#     Description of this R script:
#     R interface for multinomial sparse group lasso routines.
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

#' Compute the classification error rate 
#' 
#' @param x a msgl object
#' @param type 
#' @return a vector of error rates
#' 
#' @author Martin Vincent
#' @examples
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' 
#' # Fit model
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Training misclassification error
#' Err(fit, x)
#' 
#' # Misclassification error of xnew
#' Err(fit, x.new, classes.new)
#' 
#' # Do cross validation
#' fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Cross validation (expected) misclassification error 
#' Err(fit.cv)
#' 
#' # Do subsampling
#' test <- replicate(2, sample(1:length(classes))[1:20], simplify = FALSE)
#' train <- lapply(test, function(s) (1:length(classes))[-s])
#'
#' fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda,
#'  training = train, test = test)
#' 
#' # Mean misclassification error of the tests
#' Err(fit.sub)
#'  
#' @export
Err <- function(object, x = NULL, classes = object$classes.true, type = "rate") { #type = count, rate
	
	if(!is.null(x)) {
		return(Err(object = predict(object, x), x = NULL, type = type, classes = classes))
	}
	
	if(any(names(object) == "classes")) {
		
		if(is.null(classes)) {
			stop("true classes not found")
		}
		
		if(is.list(object$classes)) {
			#FIXME
			err.count <- lapply(object$classes, function(y) colSums(y != classes))
		} else {
			err.count <- colSums(object$classes != classes)
		}
	} else {
		stop("no classes found")				
	}
	
	if(type=="rate") {
		return(err.count/length(classes))
	}
	
	if(type=="count") {
		return(err.count)
	}
	
	stop("Unknown type")
	
}

#' todo
#' @param x 
#' @param ... 
#' @return a list of nonzero features (that is nonzero colums of the beta matrices)
#' 
#' @author martin
#' @method features msgl
#' @S3method features msgl
#' @export
features.msgl <- function(x, ...) {
	class(x) <- "sgl" # Use std function
	return(features(x))
}

#' todo
#' @param x 
#' @param ... 
#' @return todo
#' 
#' @author martin
#' @method parameters msgl
#' @S3method parameters msgl
#' @export
parameters.msgl <- function(x, ...) {
	class(x) <- "sgl" # Use std function
	return(parameters(x))
}

#' todo
#' @param x 
#' @param ... 
#' @return todo
#' 
#' @author martin
#' @method nmod msgl
#' @S3method nmod msgl
#' @export
nmod.msgl <- function(x, ...) {
	class(x) <- "sgl" # Use std function
	return(nmod(x))
}

#' Exstract the estimated models
#' 
#' @param x a msgl object 
#' @param index the models to be returned
#' @return a list of sparse matrices
#' 
#' @author Martin Vincent
#' @method models msgl
#' @S3method models msgl
#' @export
models.msgl <- function(x, index = 1:nmod(x), ...) {
	class(x) <- "sgl" # Use std function
	return(models(x))
}

#' todo
#' @param x 
#' @param index 
#' @param ... 
#' @return todo
#' 
#' @author martin
#' @method coef msgl
#' @S3method coef msgl
#' @export
coef.msgl <- function(x, index = 1:nmod(x), ...) {
	class(x) <- "sgl" # Use std function
	return(coef(x))
}


#' Print function for msgl
#'
#' This funtion will print some general information about the msgl object
#'  
#' @param x msgl object
#' @param ... not used
#' 
#' @method print msgl
#' @S3method print msgl
#' @author Martin Vincent
#' @export
print.msgl <- function(x, ...) {

	#FIXME check that  object$msgl_version exsists

        message(paste("High dimensional Multinomial logistic regression models (estimated by msgl version ", x$msgl_version, ")", sep=""))
	message()
        message(paste("  This object contains ", length(models(x)), " sparse models with ", nrow(x$beta[[1]]), " classes.", sep=""))
        message(paste("   The models contain between ", min(features(x)), " - ", max(features(x)), " nonzero features and ", min(parameters(x)), " - ", max(parameters(x)), " nonzero parameters.", sep =""))
	message()
	classnames <- paste(rownames(x$beta[[1]])[1:3], collapse=", ") #FIXME what if less than 3 classes + handle long class names
	message(paste("  The first three classes are: ", classnames, sep=""))

}

