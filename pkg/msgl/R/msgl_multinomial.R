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

#' Fit a multinomial sparse group lasso regularization path. 
#'
#' For a classification problem with  \eqn{K} classes and \eqn{p} covariates dived into \eqn{m} groups.
#' A sequence of minimizers (one for each lambda given in the \code{lambda} argument) of 
#' \deqn{\hat R(\beta) + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\hat R} is the weighted empirical log-likelihood risk of the multinomial regression model.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of covariates
#' (default is one covariate per group, hence the default dimension of \eqn{\beta^{(J)}} is \eqn{K}). 
#' The group weights \eqn{\gamma \in [0,\infty)^m} and the parameter weights \eqn{\xi = (\xi^{(1)},\dots, \xi^{(m)}) \in [0,\infty)^n} 
#' with \eqn{\xi^{(1)}\in [0,\infty)^{n_1},\dots, \xi^{(m)} \in [0,\infty)^{n_m}}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path.
#' @param return the indices of lambda values for which to return a the fitted parameters.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{beta}{the fitted parameters -- a list of length \code{length(lambda)} with each entry a matrix of size \eqn{K\times (p+1)} holding the fitted parameters}
#' \item{loss}{the values of the loss function}
#' \item{objective}{the values of the objective function (i.e. loss + penalty)}
#' \item{lambda}{the lambda values used}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#' fit$beta[[10]] #model with lambda = lambda[10] 
#' @author Martin Vincent
#' @export
#' @useDynLib msgl .registration=TRUE
#' @import Matrix
msgl <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE, 
		lambda, return = 1:length(lambda), sparse.data = FALSE, algorithm.config = sgl.standard.config) {

	# Default values
	if(is.null(grouping)) {
		covariateGrouping <- factor(1:ncol(x))
	} else {
		# ensure factor
		covariateGrouping <- factor(grouping)		
	}
	
	# cast
	classes <- factor(classes)
	return <- as.integer(return)
	
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(length(levels(classes))*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = length(levels(classes)), ncol = ncol(x))
	}
	
	# Standardize
	if(standardize) {
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- if(sparse.data) rep(0, length(x.scale)) else attr(x, "scaled:center")
	}
	
	# add intercept
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(classes))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	
	# create data
	data <- create.sgldata(x, y = NULL, sampleWeights, classes)
	
	# call SglOptimizer function
	if(data$sparseX) {
		res <- sgl_fit("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)
	} else {
		res <- sgl_fit("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, return = 1:length(lambda), algorithm.config)
	}
	
	# Convert beta back to the org scale
	if(standardize) {
		res$beta <- .to_org_scale(beta = res$beta, x.scale = x.scale, x.center = x.center)
	}
	
	class(res) <- "msgl"
	return(res)
}

#' Computes a lambda sequence for the regularization path
#' 
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda.min the smallest lambda value in the computed sequence. 
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return a vector of length \code{d} containing the compute lambda sequence.
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl .registration=TRUE
msgl.lambda.seq <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, d = 100L, standardize = TRUE, lambda.min, sparse.data = FALSE, algorithm.config = sgl.standard.config) {

	# cast
	classes <- factor(classes)
	d <- as.integer(d)
	
	# Default values
	if(is.null(grouping)) {
		covariateGrouping <- factor(1:ncol(x))
	} else {
		# ensure factor
		covariateGrouping <- factor(grouping)		
	}
		
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(length(levels(classes))*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = length(levels(classes)), ncol = ncol(x))
	}

	# Standardize
	if(standardize) {
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- if(sparse.data) rep(0, length(x.scale)) else attr(x, "scaled:center")
	}
	
	# add intercept
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(classes))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	
	# create data
	data <- create.sgldata(x, y = NULL, sampleWeights, classes)
	
	# call SglOptimizer function
	if(data$sparseX) {
		lambda <- sgl_lambda_sequence("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
	} else {
		lambda <- sgl_lambda_sequence("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha = alpha, d = d, lambda.min, algorithm.config)
	}
	
	return(lambda)
}

#' Multinomial sparse group lasso cross validation using multiple possessors 
#' 
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than \eqn{1} and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' If \code{fold}\eqn{\le}\code{max(table(classes))} then the data will be split into \code{fold} disjoint subsets keeping the ration of classes approximately equal.
#' Otherwise the data will be split into \code{fold} disjoint subsets without keeping the ration fixed.
#' @param cv.indices a list of indices of a cross validation splitting. 
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param max.threads the maximal number of threads to be used
#' @param seed the seed used for generating the random cross validation splitting, only used if \code{fold}\eqn{\le}\code{max(table(classes))}. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{link}{the linear predictors -- a list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the linear predictors.}
#' \item{response}{the estimated probabilities - a list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the probabilities.}
#' \item{classes}{the estimated classes - a matrix of size \eqn{N \times d} with \eqn{d=}\code{length(lambda)}.}
#' \item{cv.indices}{the cross validation splitting used.}
#' \item{features}{average number of features used in the models.}
#' \item{parameters}{average number of parameters used in the models.}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.03)
#' fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Missclassification count
#' colSums(fit.cv$classes != classes)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl .registration=TRUE
msgl.cv <- function(x, classes, sampleWeights = NULL, grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE, 
		lambda, fold = 10L, cv.indices = list(), sparse.data = FALSE, max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
	# Default values
	if(is.null(grouping)) {
		covariateGrouping <- factor(1:ncol(x))
	} else {
		# ensure factor
		covariateGrouping <- factor(grouping)		
	}
	
	if(is.null(sampleWeights)) {
		if(length(cv.indices) == 0) {
			sampleWeights <- rep(fold/(length(classes)*(fold-1)), length(classes))
		} else {
			n_train <- sapply(cv.indices, function(x) length(classes)-length(x))
			sampleWeights <- rep(1/mean(n_train), length(classes))
		}
	}
	
	# cast
	classes <- factor(classes)
	fold <- as.integer(fold)
	max.threads <- as.integer(max.threads)
	seed <- as.integer(seed)
		
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(length(levels(classes))*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = length(levels(classes)), ncol = ncol(x))
	}
	
	# Standardize
	if(standardize) {
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- if(sparse.data) rep(0, length(x.scale)) else attr(x, "scaled:center")
	}
	
	# add intercept
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(classes))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	
	# create data
	data <- create.sgldata(x, y = NULL, sampleWeights, classes)
		
	# call sglOptim function
	if(data$sparseX) {
		
		res <- sgl_cv("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, seed, algorithm.config)
		
	} else {
		
		res <- sgl_cv("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, seed, algorithm.config)
		
	}
	
	### Set correct dim names
	dim.names <- list(data$group.names, data$sample.names)
	
	# classes
	rownames(res$classes) <- dim.names[[2]]
	
	if(!is.null(dim.names[[1]])) {
		res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
	}
		
	# Set dim names for link and response
	res$link <- lapply(X = res$link, FUN = function(m) {dimnames(m) <- dim.names; m})
	res$response <- lapply(X = res$response, FUN = function(m) {dimnames(m) <- dim.names; m})
	
	class(res) <- "msgl"
	return(res)
}

#' Multinomial sparse group lasso generic subsampling procedure
#'
#' Support the use of multiple processors.
#' 
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path.
#' @param training a list of training samples, each item of the list corresponding to a subsample.
#' Each item in the list must be a vector with the indices of the training samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{test} list.  
#' @param test a list of test samples, each item of the list corresponding to a subsample.
#' Each item in the list must be vector with the indices of the test samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{training} list.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param max.threads the maximal number of threads to be used
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{link}{the linear predictors -- a list of length \code{length(test)} with each element of the list another list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the linear predictors.}
#' \item{response}{the estimated probabilities -- a list of length \code{length(test)} with each element of the list another list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the probabilities.}
#' \item{classes}{the estimated classes -- a list of length \code{length(test)} with each element of the list a matrix of size \eqn{N \times d} with \eqn{d=}\code{length(lambda)}.}
#' \item{features}{number of features used in the models.}
#' \item{parameters}{number of parameters used in the models.}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.03)
#'
#' test <- replicate(5, sample(1:length(classes))[1:20], simplify = FALSE)
#' train <- lapply(test, function(s) (1:length(classes))[-s])
#' 
#' fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, 
#'  training = train, test = test)
#' 
#' # Missclassification count of second subsample
#' colSums(fit.sub$classes[[2]] != classes[test[[2]]])
#' @author Martin Vincent
#' @export
#' @useDynLib msgl .registration=TRUE
msgl.subsampling <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE, 
		lambda, training, test, sparse.data = FALSE, max.threads = 2L, algorithm.config = sgl.standard.config) {
	
	# Default values
	if(is.null(grouping)) {
		covariateGrouping <- factor(1:ncol(x))
	} else {
		# ensure factor
		covariateGrouping <- factor(grouping)		
	}
	
	if(is.null(sampleWeights)) {
		if(length(cv.indices) == 0) {
			sampleWeights <- rep(fold/(length(classes)*(fold-1)), length(classes))
		} else {
			n_train <- sapply(cv.indices, function(x) length(classes)-length(x))
			sampleWeights <- rep(1/mean(n_train), length(classes))
		}
	}
	
	# cast
	classes <- factor(classes)
	max.threads <- as.integer(max.threads)
	
	if(is.null(groupWeights)) {
		groupWeights <- c(sqrt(length(levels(classes))*table(covariateGrouping)))
	}
	
	if(is.null(parameterWeights)) {
		parameterWeights <-  matrix(1, nrow = length(levels(classes)), ncol = ncol(x))
	}
	
	# Standardize
	if(standardize) {
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- if(sparse.data) rep(0, length(x.scale)) else attr(x, "scaled:center")
	}
	
	# add intercept
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	groupWeights <- c(0, groupWeights)
	parameterWeights <- cbind(rep(0, length(levels(classes))), parameterWeights)
	covariateGrouping <- factor(c("Intercept", as.character(covariateGrouping)), levels = c("Intercept", levels(covariateGrouping)))
	
	# create data
	data <- create.sgldata(x, y = NULL, sampleWeights, classes)
	
	# call sglOptim function
	if(data$sparseX) {
		
		res <- sgl_subsampling("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, max.threads, algorithm.config)
		
	} else {
		
		res <- sgl_subsampling("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, max.threads, algorithm.config)
		
	}
	
	### Reorganize

	res_reorg <- list()
	res_reorg$classes <- lapply(res$responses, function(x) x$classes + 1)
	res_reorg$response <- lapply(res$responses, function(x) x$response)
	res_reorg$link <- lapply(res$responses, function(x) x$link)
	res_reorg$features <- res$features
	res_reorg$parameters <- res$parameters
	
	res <- res_reorg
			
	### Set correct dim names
	dim.names <- list(data$group.names, data$sample.names)
	
	for(i in 1:length(test)) {
		
		#Set class names
		rownames(res$classes[[i]]) <- dim.names[[2]][test[[i]]]
		
		if(!is.null(dim.names[[1]])) {
			res$classes[[i]] <- apply(X = res$classes[[i]], MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x])
		}
		
		res$link[[i]] <- lapply(X = res$link[[i]], FUN = function(m) {dimnames(m) <- list(dim.names[[1]], dim.names[[2]][test[[i]]]); m})
		res$response[[i]] <- lapply(X = res$response[[i]], FUN = function(m) {dimnames(m) <- list(dim.names[[1]], dim.names[[2]][test[[i]]]); m})
	}
	
	class(res) <- "msgl"
	return(res)
}

.to_org_scale <- function(beta, x.scale, x.center) {
	for(l in 1:length(beta)) {
		
		beta.org <- t(t(beta[[l]])*c(1,1/x.scale))
		beta.org[,1] <- beta.org[,1] - rowSums(t(t(beta[[l]][,-1])*(x.center/x.scale)))
		
		beta[[l]] <- beta.org
	}
	
	return(beta)
}

#TODO do we need this ??
#.to_std_scale <- function(beta, x.scale, x.center) {
#	
#	for(l in 1:length(beta)) {
#		beta.std <- t(t(beta[[l]])*c(1, x.scale))
#		beta.std[,1] <- beta.std[,1] + rowSums(t(t(beta[[l]][,-1])*(x.center)))
#		
#		beta[[l]] <- beta.std
#	}
#	
#	return(beta)
#}

#' Predict
#' 
#' Computes the linear predictors, the estimated probabilities and the estimated classes for a new data set.
#'
#' @param object an object of class msgl, produced with \code{msgl}.
#' @param x a data matrix of size \eqn{N_\textrm{new} \times p}.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param ... ignored.
#' @return 
#' \item{link}{the linear predictors -- a list of length \code{length(fit$beta)} one item for each model, with each item a matrix of size \eqn{K \times N_\textrm{new}} containing the linear predictors.}
#' \item{response}{the estimated probabilities -- a list of length \code{length(fit$beta)} one item for each model, with each item a matrix of size \eqn{K \times N_\textrm{new}} containing the probabilities.}
#' \item{classes}{the estimated classes -- a matrix of size \eqn{N_\textrm{new} \times d} with \eqn{d=}\code{length(fit$beta)}.}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 20L, lambda.min = 0.01)
#' fit <- msgl(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Training error
#' res <- predict(fit, x)
#' colSums(res$classes != classes)
#' @author Martin Vincent
#' @method predict msgl
#' @S3method predict msgl
#' @export
#' @useDynLib msgl .registration=TRUE
predict.msgl <- function(object, x, sparse.data = FALSE, ...) {
	
	x <- cBind(Intercept = rep(1, nrow(x)), x)
	
	data <- list()
	
	if(is(x, "sparseMatrix")) {
		
		x <- as(x, "CsparseMatrix")
		data$X <- list(dim(x), x@p, x@i, x@x)
		
		res <- sgl_predict("msgl_sparse", "msgl", object, data)
		
	} else {
		
		data$X <- as.matrix(x)
		
		res <- sgl_predict("msgl_dense", "msgl", object, data)
		
	}
	
	### Set correct dim names
	dim.names <-  list(rownames(object$beta[[1]]), rownames(x))
	
	# classes
	rownames(res$classes) <- dim.names[[2]]
	
	if(!is.null(dim.names[[1]])) {
		res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
	}
	
	# Set dim names for link and response
	res$link <- lapply(X = res$link, FUN = function(m) {dimnames(m) <- dim.names; m})
	res$response <- lapply(X = res$response, FUN = function(m) {dimnames(m) <- dim.names; m})
	
	class(res) <- "msgl"
	return(res)
}

#' Simulated data set
#'
#' The use of this data set is only intended for testing and examples.
#' The data set contains 100 simulated samples grouped into 10 classes.
#' For each sample 400 covariates have been simulated.
#'
#' @name sim.data
#' @docType data
#' @keywords data
NULL

