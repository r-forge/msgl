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
#' @param data a matrix of 
#' @param response a vector of classes
#' @param classes a vector of classes
#' @param type type of error rate \code{rate} or \code{count}
#' @param ... not used
#' @return a vector of error rates
#' 
#' @author Martin Vincent
#' @examples
#' data(SimData)
#' x.all <- sim.data$x
#' x.1 <- sim.data$x[1:50,]
#' x.2 <- sim.data$x[51:100,]
#' classes.all <- sim.data$classes
#' classes.1 <- sim.data$classes[1:50]
#' classes.2 <- sim.data$classes[51:100]
#' 
#' #### Fit models using x.1
#' lambda <- msgl.lambda.seq(x.1, classes.1, alpha = .5, d = 25L, lambda.min = 0.03)
#' fit <- msgl(x.1, classes.1, alpha = .5, lambda = lambda)
#' 
#' #### Training errors:
#' 
#' # Misclassification rate
#' Err(fit, x.1)
#' 
#' # Misclassification count
#' Err(fit, x.1, type = "count")
#' 
#' # Negative log likelihood error
#' Err(fit, x.1, type="loglike")
#'  
#' # Misclassification rate of x.2
#' Err(fit, x.2, classes.2)
#' 
#' #### Do cross validation
#' fit.cv <- msgl.cv(x.all, classes.all, alpha = .5, lambda = lambda)
#' 
#' #### Cross validation errors (estimated expected generalization error)
#' 
#' # Misclassification rate
#' Err(fit.cv)
#' 
#' # Negative log likelihood error
#' Err(fit.cv, type="loglike")
#' 
#' #### Do subsampling
#' test <- list(1:20, 21:40)
#' train <- lapply(test, function(s) (1:length(classes.all))[-s])
#'
#' fit.sub <- msgl.subsampling(x.all, classes.all, alpha = .5, lambda = lambda, training = train, test = test)
#' 
#' # Mean misclassification error of the tests
#' Err(fit.sub)
#' 
#' # Negative log likelihood error
#' Err(fit.sub, type="loglike")
#'  
#' @method Err msgl
#' @S3method Err msgl
#' @import sglOptim
#' @export
Err.msgl <- function(x, data = NULL, response = x$classes.true, classes = response, type = "rate", ... ) {
	
	if(type=="rate") {
		return(compute_error(x, data = data, response.name = "classes", response = response, loss = function(x,y) mean(x != y)))
	}
	
	if(type=="count") {
		return(compute_error(x, data = data, response.name = "classes", response = response, loss = function(x,y) sum(x != y)))
	}
	
	if(type=="loglike") {
		loss <- function(x,y) -mean(log(sapply(1:length(y), function(i) x[as.integer(y[i]),i])))
		return(compute_error(x, data = data, response.name = "response", response = response, loss = loss))
	}
	
	stop("Unknown type")
	
}

#' todo
#' @param x a msgl object
#' @param ... not used
#' @return a list of nonzero features (that is nonzero colums of the beta matrices)
#' 
#' @author martin
#' @method features msgl
#' @S3method features msgl
#' @import sglOptim
#' @export
features.msgl <- function(x, ...) {
	class(x) <- "sgl" # Use std function
	return(features(x))
}

#' todo
#' @param x a msgl object
#' @param ... not used
#' @return todo
#' 
#' @author martin
#' @method parameters msgl
#' @S3method parameters msgl
#' @import sglOptim
#' @export
parameters.msgl <- function(x, ...) {
	class(x) <- "sgl" # Use std function
	return(parameters(x))
}

#' todo
#' @param x a msgl object
#' @param ... not used
#' @return todo
#' 
#' @author martin
#' @method nmod msgl
#' @S3method nmod msgl
#' @import sglOptim
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
#' @import sglOptim
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
#' @import sglOptim
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
#' @import sglOptim
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
