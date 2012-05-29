# TODO: Add comment
# 
# Author: martin
###############################################################################

#' Computes a lambda sequence for the logistic sparse group lasso regularization path
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param group.dim vector of group dimensions 
#' @param groupWeights group weights
#' @param parameterWeights parameter weights
#' @param alpha 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty
#' @param d length of the required sequence 
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda.min minimum lambda value
#' @param sparse.data if TRUE x will be treated as sparse
#' @param algorithm.config the algorithm configuration to be used 
#' @returnType 
#' @return lambda sequence
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes > 5
#' 
#' group.dim <- as.integer(c(1, rep(10, 40)))
#' group.dim <- as.integer(group.dim)
#' group.weights <- rep(1, length(group.dim))
#' parameter.weights <- rep(1, 401)
#' group.weights[1] <- 0 # no penalty on the intercept
#' parameter.weights[1] <- 0 
#' 
#' lambda <- logreg.sgl.lambda.seq(x, classes, group.dim, group.weights, parameter.weights, lambda.min = 0.01)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_logreg_sgl_lambda_seq
logreg.sgl.lambda.seq <- function(x, classes, group.dim, groupWeights, parameterWeights, alpha = 0.5, d = 100L, standardize = TRUE, lambda.min, sparse.data = FALSE, algorithm.config = sgl.standard.config) {
	
	#TODO factor(classes) will not work on a numeric vector - 
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}

	if(is.list(classes) && length(classes) == 1) {
		classes <- classes[[1]]
	}
	
	if(is.list(classes)) {
		classes.numeric <- lapply(classes, function(x) as.integer(factor(x)) - 1L)
		classes.numeric <- lapply(classes.numeric, function(x) { x[is.na(x)] <- -1L; x })
	} else if(is.vector(classes) || is.factor(classes)) {
		classes.numeric <- as.integer(factor(classes))-1L
	} else {
		stop("classes is of unsupported type.")
	}
	
	if(standardize) {
		
		if(sparse.data) {
			x <- scale(x, FALSE, TRUE)
			
		} else {	
			x <- scale(x, TRUE, TRUE)
			
		}
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	if(is.list(classes)) {
		
		stop("Support for tensor mode not implemented")
		
	} else  {
		
		if(sparse.data) {
	
			stop("Sparse matrix support not yet implemented")
				
		} else {
			res <- .Call(r_logreg_sgl_lambda_seq, x, classes.numeric, group.dim, groupWeights, parameterWeights, alpha, d, lambda.min, algorithm.config)
		}
	}
	
	return(res)
}

#' Fit a logistic sparse group lasso regularization path
#' 
#' Note that, if standardize = TRUE the fitted parameters are currently not rescaled to orgnial scale before returned.
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param group.dim vector of group dimensions 
#' @param groupWeights group weights
#' @param parameterWeights parameter weights
#' @param alpha 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path
#' @param return the indices of lambda values for which to return a the fitted parameters
#' @param do.refit if TRUE a refitted model will be returned
#' @param sparse.data if TRUE x will be treated as sparse
#' @param algorithm.config the algorithm configuration to be used 
#' @returnType 
#' @return 
#' \item{beta}{the fittede paramters}
#' \item{loss}{the values of the loss function}
#' \item{objective}{the values of the objective function (i.e. loss + penalty)}
#' \item{lambda}{the lambda values used}
#' \item{beta.refit}{refitted paramters (only if do_refit = TRUE)}
#' \item{loss.refit}{the values of the refitted loss function (only if do_refit = TRUE)}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes > 5
#' 
#' group.dim <- as.integer(c(1, rep(10, 40)))
#' group.dim <- as.integer(group.dim)
#' group.weights <- rep(1, length(group.dim))
#' parameter.weights <- rep(1, 401)
#' group.weights[1] <- 0 # no penalty on the intercept
#' parameter.weights[1] <- 0 
#' 
#' lambda <- logreg.sgl.lambda.seq(x, classes, group.dim, group.weights, parameter.weights, lambda.min = 0.01)
#' fit <- logreg.sgl(x, classes, group.dim, group.weights, parameter.weights, lambda = lambda, standardize = FALSE)
#' fit$beta[[10]] #model with lambda = lambda[10] 
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_logreg_sgl_basic
logreg.sgl <- function(x, classes, group.dim, groupWeights, parameterWeights, alpha = 0.5, standardize = TRUE, 
		lambda, return = 1:length(lambda), do.refit = FALSE, sparse.data = FALSE, algorithm.config = sgl.standard.config) {
	
	#TODO check that return is valid
	return <- as.integer(sort(unique(return))) - 1L
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	if(is.list(classes) && length(classes) == 1) {
		classes <- classes[[1]]
	}
	
	if(is.list(classes)) {
		classes.numeric <- lapply(classes, function(x) as.integer(factor(x)) - 1L)
		classes.numeric <- lapply(classes.numeric, function(x) { x[is.na(x)] <- -1L; x })
	} else if(is.vector(classes) || is.factor(classes)) {
		classes.numeric <- as.integer(factor(classes))-1L
	} else {
		stop("classes is of unsupported type.")
	}
	
	if(standardize) {
		
		if(sparse.data) {
			
			x <- scale(x, FALSE, TRUE)
			x.scale <- attr(x,"scaled:scale")
			x.center <- rep(0, length(x.scale))
			
		} else {
			
			x <- scale(x, TRUE, TRUE)
			x.scale <- attr(x, "scaled:scale")
			x.center <- attr(x, "scaled:center")
		}
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	# Run msgl algorithm
	if(is.list(classes)) {
		
		stop("Support for tensor mode not implemented")
				
	} else {
		
		if(sparse.data) {
			
			stop("Sparse matrix support not yet implemented")
			
		} else {
			
			res <- .Call(r_logreg_sgl_basic, x, classes.numeric, group.dim, groupWeights, parameterWeights, alpha, lambda, return, do.refit, algorithm.config);
			
		}
	}
	
	#TODO rescale beta
#	# Dim names
#	feature.names <- if(!is.null(colnames(x))) colnames(x) else 1:dim(x)[2]
#	
#	if(is.list(classes)) {
#		
#		class.names <- unlist(lapply(classes, function(x) levels(factor(x))))
#		
#	} else {
#		
#		class.names <- levels(factor(classes))
#		
#	}
#	
#	# Create R sparse matrix
#	res$beta <- lapply(1:length(res$beta), function(i) sparseMatrix(i = res$beta[[i]][[2]], j = res$beta[[i]][[3]], x = res$beta[[i]][[4]], dims = res$beta[[i]][[1]], dimnames = list(class.names, c("Intercept", feature.names)), index1 = FALSE))
#	
#	if(do.refit) {
#		res$beta.refit <- lapply(1:length(res$beta.refit), function(i) sparseMatrix(i = res$beta.refit[[i]][[2]], j = res$beta.refit[[i]][[3]], x = res$beta.refit[[i]][[4]], dims = res$beta.refit[[i]][[1]], dimnames = list(class.names, c("Intercept", feature.names)), index1 = FALSE))
#	}
#	
#	#Dimnames on cirtical.bounds
#	res$critical.bounds <- lapply(res$critical.bounds, function(x) {names(x) <- c("Intercept", feature.names); return(x)})
#	
#	# Convert beta back to the org scale
#	if(standardize) {
#		
#		res$beta <- .to_org_scale(beta = res$beta, x.scale = x.scale, x.center = x.center)
#		
#		if(do.refit) {
#			res$beta.refit <- .to_org_scale(beta = res$beta.refit, x.scale = x.scale, x.center = x.center)
#		}
#		
#	}
#	
#	#TODO name
#	res$classes.prior <- classes

	class(res) <- "lrsgl"
	return(res)
}

#' Logistic sparse group lasso cross validation
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param group.dim vector of group dimensions 
#' @param groupWeights group weights
#' @param parameterWeights parameter weights
#' @param alpha 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path
#' @param fold 
#' @param cv.indices 
#' @param do.refit if TRUE a refitted model will be returned
#' @param sparse.data if TRUE x will be treated as sparse
#' @param max.threads maximal number of threads
#' @param seed 
#' @param algorithm.config the algorithm configuration to be used 
#' @returnType 
#' @return 
#' \item{link}{linear predictors}
#' \item{response}{estimated probabilities}
#' \item{classes}{estimated classes}
#' \item{link.refit}{linear predictors for the refitted models (only if do_refit = TRUE)}
#' \item{response.refit}{estimated probabilities for the refitted models (only if do_refit = TRUE)}
#' \item{classes.refit}{estimated classes for the refitted models (only if do_refit = TRUE)}
#' \item{cv.indices}{}
#' \item{groups}{Average number of groups used in the models}
#' \item{parameters}{Average number of parameters used in the models}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes > 5
#' 
#' group.dim <- as.integer(c(1, rep(10, 40)))
#' group.dim <- as.integer(group.dim)
#' group.weights <- rep(1, length(group.dim))
#' parameter.weights <- rep(1, 401)
#' group.weights[1] <- 0 # no penalty on the intercept
#' parameter.weights[1] <- 0 
#' 
#' lambda <- logreg.sgl.lambda.seq(x, classes, group.dim, group.weights, parameter.weights, lambda.min = 0.01)
#' fit.cv <- logreg.sgl.cv(x, classes, group.dim, group.weights, parameter.weights, lambda = lambda, fold = 10L, max.threads = 2L)
#' 
#' # Missclassification count
#' colSums(fit.cv$classes != classes)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_logreg_sgl_cv
logreg.sgl.cv <- function(x, classes, group.dim, groupWeights, parameterWeights, alpha = 0.5, standardize = TRUE, 
		lambda, fold = 10L, cv.indices = list(), do.refit = FALSE, sparse.data = FALSE, max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	if(is.list(classes) && length(classes) == 1) {
		classes <- classes[[1]]
	}
	
	if(is.list(classes)) {
		classes.numeric <- lapply(classes, function(x) as.integer(factor(x)) - 1L)
		classes.numeric <- lapply(classes.numeric, function(x) { x[is.na(x)] <- -1L; x })
	} else if(is.vector(classes) || is.factor(classes)) {
		classes.numeric <- as.integer(factor(classes))-1L
	} else {
		stop("classes is of unsupported type.")
	}
	
	if(standardize) {
		
		if(sparse.data) {
			x <- scale(x, FALSE, TRUE)
			
		} else {	
			x <- scale(x, TRUE, TRUE)
			
		}
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	#Dimnames
	
	if(is.list(classes)) {
		class.names <- unlist(lapply(classes, function(x) levels(factor(x))))
	} else {
		class.names <- levels(factor(classes))
	}
	
	dim.names <-  list(class.names, rownames(x))
	
	if(is.list(classes)) {
		
		stop("Not implemented")
		
		
	} else {
		
		if(sparse.data) {
			
			stop("Not implemented")
			
		} else {
			
			if(length(cv.indices) == 0) {
				res <- .Call(r_logreg_sgl_cv, x, classes.numeric, group.dim, groupWeights, parameterWeights, alpha, lambda, do.refit, fold, cv.indices, FALSE, max.threads, seed, algorithm.config)
			} else {
				cv.indices <- lapply(cv.indices, function(x) as.integer(x-1))
				res <- .Call(r_logreg_sgl_cv, x, classes.numeric, group.dim, groupWeights, parameterWeights, alpha, lambda, do.refit, fold, cv.indices, TRUE, max.threads, seed, algorithm.config)
			}
		}
		
		#Set class names
		rownames(res$classes) <- dim.names[[2]]
		
		if(!is.null(dim.names[[1]])) {
			res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
		}
		
		if(do.refit) {
			rownames(res$classes.refit) <- dim.names[[2]]
			
			if(!is.null(dim.names[[1]])) {
				res$classes.refit <- apply(X = res$classes.refit, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
			}
		}
		
	}
	
	res$link <- lapply(X = res$link, FUN = function(m) {dimnames(m) <-dim.names; m})
	res$response <- lapply(X = res$response, FUN = function(m) {dimnames(m) <-dim.names; m})
	
	class(res) <- "lrsgl"
	return(res)
}

#' Predict
#' 
#' @param fit object of class 'lrsgl'
#' @param x new data
#' @param sparse.data if TRUE x will be treated as sparse
#' @returnType 
#' @return 
#' \item{link}{linear predictors}
#' \item{response}{estimated probabilities}
#' \item{classes}{estimated classes}
#' @author Martin Vincent
#' @export
#' @useDynLib msgl
predict.lrsgl <- function(fit, x, sparse.data = FALSE) {
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	res <- .predict.lrsgl(fit$beta, x, sparse.data, is.list(fit$classes))	
	
	if("beta.refit" %in% names(fit)) {
		res.refit <- .predict.lrsgl(fit$beta.refit, x)
		
		res$classes.refit <- res.refit$classes
		res$response.refit <- res.refit$response
		res$link.refit <- res.refit$link
	}
	
	return(res)
}

.predict.lrsgl <- function(beta, x, sparse.data = FALSE, tensor = FALSE) {
	
	#TODO dimnames

	classes <- fit$classes.prior
	
	if(!tensor) {
		
		if(sparse.data) {
			
			stop("Not implemented")
			
		} else {
			
			res <- .Call(r_logreg_predict, x, beta)
		}

	} else {

		stop("Not implemented")

	}
	
	return(res)
}

