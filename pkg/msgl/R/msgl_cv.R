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
                lambda, fold = 10L, cv.indices = list(), sparse.data = FALSE, max.threads = 2L, algorithm.config = sgl.standard.config) {

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

                res <- sgl_cv("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, algorithm.config)

        } else {

                res <- sgl_cv("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, max.threads, algorithm.config)

        }

		### Responses
		res$classes <- res$responses$classes
		res$response <- res$responses$response
		res$link <- res$responses$link
		res$responses <- NULL
		
       # Set class names
		if(!is.null(data$group.names)) {
			res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) data$group.names[x])
			res$link <- lapply(X = res$link, FUN = function(m) {rownames(m) <- data$group.names; m})
			res$response <- lapply(X = res$response, FUN = function(m) {rownames(m) <- data$group.names; m})
		}

        res$msgl_version <- packageVersion("msgl")

        class(res) <- "msgl"
        return(res)
}