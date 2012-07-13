# TODO: Add comment
# 
# Author: martin
###############################################################################

#TODO classWeights when classes is list

#' Fit multinomial sparse group lasso regularization path
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param featureWeights the group weights, a vector of length p
#' @param classWeights a vector of length K
#' @param alpha 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path
#' @param return the indices of lambda values for which to return a the fitted parameters
#' @param do.refit if TRUE a refitted model will be returned
#' @param sparse.data if TRUE x will be treated as sparse
#' @param algorithm.config the algorithm configuration to be used 
#' @returnType S3 object of class 'msgl' 
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
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' fit <- msgl(x, classes, lambda = lambda)
#' fit$beta[[10]] #model with lambda = lambda[10] 
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_sparse_basic r_msgl_basic
msgl <- function(x, classes, featureWeights = .featureWeights(x, classes), classWeights = .classWeights(classes), alpha = 0.5, standardize = TRUE, 
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
		
		if(sparse.data) {
			stop("Sparse matrix support not yet implemented for tensor mode")
		}
		
		res <- .Call(r_msgl_tensor_basic, x, classes.numeric, featureWeights, classWeights, alpha, lambda, return, do.refit, algorithm.config);
		
	} else {
		
		if(sparse.data) {
			
			m <- as(x, "TsparseMatrix")
			m <- list(dim(m), m@i, m@j, m@x)
			
			res <- .Call(r_msgl_sparse_basic, m, classes.numeric, featureWeights, classWeights, alpha, lambda, return, do.refit, algorithm.config);
			
		} else {
			
			res <- .Call(r_msgl_basic, x, classes.numeric, featureWeights, classWeights, alpha, lambda, return, do.refit, algorithm.config);
			
		}
	}
	
	# Dim names
	feature.names <- if(!is.null(colnames(x))) colnames(x) else 1:dim(x)[2]
	
	if(is.list(classes)) {
		
		class.names <- unlist(lapply(classes, function(x) levels(factor(x))))
		
	} else {
		
		class.names <- levels(factor(classes))
		
	}
	
	# Create R sparse matrix
	res$beta <- lapply(1:length(res$beta), function(i) sparseMatrix(i = res$beta[[i]][[2]], j = res$beta[[i]][[3]], x = res$beta[[i]][[4]], dims = res$beta[[i]][[1]], dimnames = list(class.names, c("Intercept", feature.names)), index1 = FALSE))
	
	if(do.refit) {
		res$beta.refit <- lapply(1:length(res$beta.refit), function(i) sparseMatrix(i = res$beta.refit[[i]][[2]], j = res$beta.refit[[i]][[3]], x = res$beta.refit[[i]][[4]], dims = res$beta.refit[[i]][[1]], dimnames = list(class.names, c("Intercept", feature.names)), index1 = FALSE))
	}
	
	#Dimnames on cirtical.bounds
	res$critical.bounds <- lapply(res$critical.bounds, function(x) {names(x) <- c("Intercept", feature.names); return(x)})
	
	# Convert beta back to the org scale
	if(standardize) {
		
		res$beta <- .to_org_scale(beta = res$beta, x.scale = x.scale, x.center = x.center)
		
		if(do.refit) {
			res$beta.refit <- .to_org_scale(beta = res$beta.refit, x.scale = x.scale, x.center = x.center)
		}
		
	}
	
	#TODO name
	res$classes.prior <- classes
	
	class(res) <- "msgl"; return(res)
}

#' Computes a lambda sequence for the regularization path
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param featureWeights the group weights, a vector of length p
#' @param classWeights a vector of length K
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
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_sparse_lambda_seq r_msgl_lambda_seq
msgl.lambda.seq <- function(x, classes, featureWeights = .featureWeights(x, classes), classWeights = .classWeights(classes), alpha = 0.5, d = 100L, standardize = TRUE, lambda.min, sparse.data = FALSE, algorithm.config = sgl.standard.config) {
	
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
		
		if(sparse.data) {
			stop("Sparse matrix support not yet implemented for tensor mode")
		}
		
		res <- .Call(r_msgl_tensor_lambda_seq, x, classes.numeric, featureWeights, classWeights, alpha, d, lambda.min, algorithm.config)
	} else  {
		
		if(sparse.data) {
			
			x <- as(x, "TsparseMatrix")
			x <- list(dim(x), x@i, x@j, x@x)
			
			res <- .Call(r_msgl_sparse_lambda_seq, x, classes.numeric, featureWeights, classWeights, alpha, d, lambda.min, algorithm.config)
			
		} else {
			res <- .Call(r_msgl_lambda_seq, x, classes.numeric, featureWeights, classWeights, alpha, d, lambda.min, algorithm.config)
		}
	}
	
	return(res)
}

#' Multinomial sparse group lasso cross validation
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param featureWeights the group weights, a vector of length p
#' @param classWeights a vector of length K
#' @param alpha 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path
#' @param fold 
#' @param cv.classes 
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
#' \item{features}{Average number of features used in the models}
#' \item{parameters}{Average number of parameters used in the models}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Missclassification count
#' colSums(fit.cv$classes != classes)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_sparse_cv r_msgl_cv
msgl.cv <- function(x, classes, featureWeights = .featureWeights(x, classes), classWeights = .classWeights(classes), alpha = 0.5, standardize = TRUE, 
		lambda, fold = 10L, cv.classes = .list.factor(classes), cv.indices = list(), do.refit = FALSE, sparse.data = FALSE, max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
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
		
		if(sparse.data) {
			stop("Sparse matrix support not yet implemented for tensor mode")
		}
		
		if(length(cv.indices) == 0) {
			res <- .Call(r_msgl_tensor_cv, x, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, fold, cv.classes, cv.indices, FALSE, max.threads, seed, algorithm.config)
		} else {
			cv.indices <- lapply(cv.indices, function(x) as.integer(x-1))
			res <- .Call(r_msgl_tensor_cv, x, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, fold, cv.classes, cv.indices, TRUE, max.threads, seed, algorithm.config)
		}
		
		res$classes <- lapply(X = res$classes, FUN = function(m) {
					
					m.new <- matrix(nrow = nrow(m), ncol = ncol(m))
					
					colnames(m.new) <- dim.names[[2]]; 
					
					#TODO faster implemtation 
					for(i in 1:nrow(m)) {
						for(j in 1:ncol(m)) {
							m.new[i,j] <- levels(classes[[i]])[m[i,j]+1]
						}
					}
					
					return(m.new)
				}) 
		
		if(do.refit) {
			
			res$classes.refit <- lapply(X = res$classes.refit, FUN = function(m) {
						
						m.new <- matrix(nrow = nrow(m), ncol = ncol(m))
						
						colnames(m.new) <- dim.names[[2]]; 
						
						#TODO faster implemtation 
						for(i in 1:nrow(m)) {
							for(j in 1:ncol(m)) {
								m.new[i,j] <- levels(classes[[i]])[m[i,j]+1]
							}
						}
						
						return(m.new)
					}) 
		}
		
	} else {
		
		if(sparse.data) {
			
			m <- as(x, "TsparseMatrix")
			m <- list(dim(m), m@i, m@j, m@x)
			
			if(length(cv.indices) == 0) {
				res <- .Call(r_msgl_sparse_cv, m, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, fold, cv.indices, FALSE, max.threads, seed, algorithm.config)
			} else {
				cv.indices <- lapply(cv.indices, function(x) as.integer(x-1))
				res <- .Call(r_msgl_sparse_cv, m, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, fold, cv.indices, TRUE, max.threads, seed, algorithm.config)
			}
			
		} else {
			
			if(length(cv.indices) == 0) {
				res <- .Call(r_msgl_cv, x, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, fold, cv.indices, FALSE, max.threads, seed, algorithm.config)
			} else {
				cv.indices <- lapply(cv.indices, function(x) as.integer(x-1))
				res <- .Call(r_msgl_cv, x, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, fold, cv.indices, TRUE, max.threads, seed, algorithm.config)
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
	
	
	class(res) <- "msgl"
	return(res)
}


#'  Predict
#' 
#' @param fit object of class 'msgl'
#' @param x new data
#' @param sparse.data if TRUE x will be treated as sparse
#' @returnType 
#' @return 
#' \item{link}{linear predictors}
#' \item{response}{estimated probabilities}
#' \item{classes}{estimated classes}
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_predict_classes r_msgl_sparse_predict_classes
predict.msgl <- function(fit, x, sparse.data = FALSE) {
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	res <- .predict.msgl(fit$beta, x, sparse.data, is.list(fit$classes))	
	
	if("beta.refit" %in% names(fit)) {
		res.refit <- .predict.msgl(fit$beta.refit, x)
		
		res$classes.refit <- res.refit$classes
		res$response.refit <- res.refit$response
		res$link.refit <- res.refit$link
	}
	
	class(res) <- "msgl"
	return(res)
}

.predict.msgl <- function(beta, x, sparse.data = FALSE, tensor = FALSE) {
	#Save dimnames
	dim.names <-  list(rownames(beta[[1]]), rownames(x))
	
	beta <- lapply(X = beta, FUN = function(m) as(m, "TsparseMatrix"))
	beta <- lapply(X = beta, FUN = function(m) list(dim(m), m@i, m@j, m@x))
	
	classes <- fit$classes.prior
	
	if(!tensor) {
		
		if(sparse.data) {
			
			m <- as(x, "TsparseMatrix")
			m <- list(dim(m), m@i, m@j, m@x)
			
			res <- .Call(r_msgl_sparse_predict_classes, m, beta)
			
		} else {
					
			res <- .Call(r_msgl_predict_classes, x, beta)
		}
		
		#Set class names
		rownames(res$classes) <- dim.names[[2]]
		
		if(!is.null(dim.names[[1]])) {
			res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
		}
		
	} else {
		
		if(sparse.data) {
			stop("Sparse matrix support not yet implemented for tensor mode")
		}
		
		res <- .Call(r_msgl_tensor_predict_classes, x, beta, sapply(classes, function(x) length(levels(x))))
		
		res$classes <- lapply(X = res$classes, FUN = function(m) {
					
					m.new <- matrix(nrow = nrow(m), ncol = ncol(m))
					
					colnames(m.new) <- dim.names[[2]]; 
					
					#TODO faster implemtation 
					for(i in 1:nrow(m)) {
						for(j in 1:ncol(m)) {
							m.new[i,j] <- levels(classes[[i]])[m[i,j]+1]
						}
					}
					
					return(m.new)
				}) 
	}
	
	
	#Set dimnames ect
	
	res$link <- lapply(X = res$link, FUN = function(m) {dimnames(m) <-dim.names; m})
	res$response <- lapply(X = res$response, FUN = function(m) {dimnames(m) <-dim.names; m})
	
	
	class(res) <- "msgl"
	return(res)
}

#' Multinomial sparse group lasso subsampling
#' 
#' @param x design matrix, matrix of size N x p
#' @param classes grouping of samples, factor of length N
#' @param featureWeights the group weights, a vector of length p
#' @param classWeights a vector of length K
#' @param alpha 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path
#' @param cv.classes unused
#' @param d number of subsamples
#' @param fraction
#' @param subsamples alternative way of specifing the subsamples. A list of index vectors. 
#' @param do.refit if TRUE a refitted model will be returned
#' @param sparse.data if TRUE x will be treated as sparse
#' @param max.threads maximal number of threads
#' @param seed 
#' @param algorithm.config the algorithm configuration to be used 
#' @returnType 
#' @return ...
#' @author Martin vincent
#' @export
msgl.subsampleing <- function(x, classes, featureWeights = .featureWeights(x, classes), classWeights = .classWeights(classes), alpha = 0.5, standardize = TRUE, 
		lambda, cv.classes = .list.factor(classes), d = 100L, fraction = 0.2, subsamples = list(), do.refit = FALSE, sparse.data = FALSE, max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
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
			
			if(length(subsamples) == 0) {
				res <- .Call(r_msgl_subsampleing, x, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, d, fraction, subsamples, FALSE, max.threads, seed, algorithm.config)
			} else {
				subsamples <- lapply(subsamples, function(x) as.integer(x - 1))
				res <- .Call(r_msgl_subsampleing, x, classes.numeric, featureWeights, classWeights, alpha, lambda, do.refit, d, fraction, subsamples, TRUE, max.threads, seed, algorithm.config)
			}
		}
	}
	
	# Convert 0 based indexing to 1 based
	res$subsamples <- lapply(X = res$subsamples, FUN = function(x) as.integer(x+1))	

	#Set class names
	
	subsamples <- res$subsamples

	for(i in 1:length(subsamples)) {
	
		rownames(res$classes[[i]]) <- dim.names[[2]][subsamples[[i]]]
		
		if(!is.null(dim.names[[1]])) {
			res$classes[[i]] <- apply(X = res$classes[[i]], MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
		}
		
		if(do.refit) {
			rownames(res$classes.refit[[i]]) <- dim.names[[2]][subsamples[[i]]]
			
			if(!is.null(dim.names[[1]])) {
				res$classes.refit[[i]] <- apply(X = res$classes.refit[[i]], MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
			}
		}

		res$link[[i]] <- lapply(X = res$link[[i]], FUN = function(m) {dimnames(m) <- list(dim.names[[1]], dim.names[[2]][subsamples[[i]]]); m})
		res$response[[i]] <- lapply(X = res$response[[i]], FUN = function(m) {dimnames(m) <- list(dim.names[[1]], dim.names[[2]][subsamples[[i]]]); m})
		
	}

	class(res) <- "msgl"
	return(res)
}


#' Create a new algorithm configuration
#' 
#' @param tolerance_penalized_main_equation_loop 
#' @param tolerance_penalized_inner_loop_alpha 
#' @param tolerance_penalized_inner_loop_beta 
#' @param tolerance_penalized_middel_loop_alpha 
#' @param tolerance_penalized_middel_loop_beta 
#' @param tolerance_penalized_outer_loop_alpha 
#' @param tolerance_penalized_outer_loop_beta 
#' @param tolerance_penalized_outer_loop_gamma 
#' @param general_tolerance_unpenalized 
#' @param use_bound_optimization 
#' @param use_active_set_optimization 
#' @param active_set_threshold_alpha 
#' @param active_set_threshold_beta 
#' @param use_stepsize_optimization_in_penalizeed_loop 
#' @param stepsize_opt_penalized_initial_t 
#' @param stepsize_opt_penalized_a 
#' @param stepsize_opt_penalized_b 
#' @param stepsize_opt_unpenalized_initial_t 
#' @param stepsize_opt_unpenalized_a 
#' @param stepsize_opt_unpenalized_b 
#' @param verbose 
#' @returnType 
#' @return ...
#' @author Martin Vincent
#' @export
sgl.algorithm.config <- function(tolerance_penalized_main_equation_loop = 1e-10, tolerance_penalized_inner_loop_alpha = 1e-4, tolerance_penalized_inner_loop_beta = 5, tolerance_penalized_middel_loop_alpha = 0.01, tolerance_penalized_middel_loop_beta = 0,
		tolerance_penalized_outer_loop_alpha = 0.01, tolerance_penalized_outer_loop_beta = 0, tolerance_penalized_outer_loop_gamma = 1e-5, general_tolerance_unpenalized = 1e-3, use_bound_optimization = TRUE, 
		use_active_set_optimization = FALSE, active_set_threshold_alpha = 5e-1, active_set_threshold_beta = 0,
		use_stepsize_optimization_in_penalizeed_loop = TRUE, stepsize_opt_penalized_initial_t = 1,
		stepsize_opt_penalized_a = 0.1, stepsize_opt_penalized_b = 0.5, stepsize_opt_unpenalized_initial_t = 1,
		stepsize_opt_unpenalized_a = 0.1, stepsize_opt_unpenalized_b = 0.5, verbose = FALSE) {
	
	config <- list()
	
	config$tolerance_penalized_main_equation_loop <- tolerance_penalized_main_equation_loop
	
	config$tolerance_penalized_inner_loop_alpha <- tolerance_penalized_inner_loop_alpha
	config$tolerance_penalized_inner_loop_beta <- tolerance_penalized_inner_loop_beta
	
	config$tolerance_penalized_middel_loop_alpha <- tolerance_penalized_middel_loop_alpha
	config$tolerance_penalized_middel_loop_beta <- tolerance_penalized_middel_loop_beta
	
	config$tolerance_penalized_outer_loop_alpha <- tolerance_penalized_outer_loop_alpha
	config$tolerance_penalized_outer_loop_beta <- tolerance_penalized_outer_loop_beta	
	config$tolerance_penalized_outer_loop_gamma <- tolerance_penalized_outer_loop_gamma	
	
	config$general_tolerance_unpenalized <- general_tolerance_unpenalized
	
	config$use_bound_optimization <- use_bound_optimization
	
	config$use_active_set_optimization <- use_active_set_optimization
	config$active_set_threshold_alpha <- active_set_threshold_alpha
	config$active_set_threshold_beta <- active_set_threshold_beta
	
	config$use_stepsize_optimization_in_penalizeed_loop <- use_stepsize_optimization_in_penalizeed_loop
	config$stepsize_opt_penalized_initial_t <- stepsize_opt_penalized_initial_t
	config$stepsize_opt_penalized_a <- stepsize_opt_penalized_a
	config$stepsize_opt_penalized_b <- stepsize_opt_penalized_b
	
	config$stepsize_opt_unpenalized_initial_t <- stepsize_opt_unpenalized_initial_t
	config$stepsize_opt_unpenalized_a <- stepsize_opt_unpenalized_a
	config$stepsize_opt_unpenalized_b <- stepsize_opt_unpenalized_b
	
	config$verbose <- verbose
	
	return(config)
}

#' Standard algorithm configuartion
#' 
#' @author Martin Vicnet
#' @export
sgl.standard.config <- sgl.algorithm.config();

.to_org_scale <- function(beta, x.scale, x.center) {
	for(l in 1:length(beta)) {
		
		beta.org <- t(t(beta[[l]])*c(1,1/x.scale))
		beta.org[,1] <- beta.org[,1] - rowSums(t(t(beta[[l]][,-1])*(x.center/x.scale)))
		
		beta[[l]] <- beta.org
	}
	
	return(beta)
}

.to_std_scale <- function(beta, x.scale, x.center) {
	
	for(l in 1:length(beta)) {
		beta.std <- t(t(beta[[l]])*c(1, x.scale))
		beta.std[,1] <- beta.std[,1] + rowSums(t(t(beta[[l]][,-1])*(x.center)))
		
		beta[[l]] <- beta.std
	}
	
	return(beta)
}

.list.factor <- function(l) {
	
	if(!is.list(l)) {
		return (as.integer(l)-1L)
	}
	
	m <- sapply(l, function(x) as.character(x))
	
	v <- vector("integer", nrow(m))
	u <- unique(m)
	
	for(i in 1:nrow(u)) {
		v[v==0] = apply(m, MARGIN = 1, function(x) all(x == u[i,], na.rm=TRUE)*i)[v==0]
	}
	
	v = v - 1L;
	
	return(v)
}

.count.levels <- function(x) {
	
	if(is.list(x)) {
		return(sum(sapply(x, function(x) length(levels(x)))))
	}
	
	return(length(levels(factor(x))))
}

.featureWeights <- function(x, classes) {
	
	if(is.list(x)) { #TODO is(x, "sparseMatrix")) will not work, it seems R do not know that x is a s4 object
		return(rep(sqrt(.count.levels(classes)), x[[1]][2])) #TODO a hack
	} else {
		return(rep(sqrt(.count.levels(classes)), ncol(x)))
	}
	
}

.classWeights <- function(classes) {
	return(rep(1, .count.levels(classes)))
}
