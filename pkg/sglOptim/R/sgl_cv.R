#
#     Description of this R script:
#     R interface generic sparse group lasso cross validation
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

#' Generic sparse group lasso cross validation using multiple possessors 
#' 
#' 
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}. 
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups). 
#' @param parameterWeights a matrix of size \eqn{q \times p}. 
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than \eqn{1} and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' If \code{fold}\eqn{\le}\code{max(table(classes))} then the data will be split into \code{fold} disjoint subsets keeping the ration of classes approximately equal.
#' Otherwise the data will be split into \code{fold} disjoint subsets without keeping the ration fixed.
#' @param cv.indices a list of indices of a cross validation splitting. 
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param max.threads the maximal number of threads to be used.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return sgl object content will depend on the C++ response class.
#' @author Martin Vincent
#' @export
sgl_cv <- function(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, lambda, fold = 2L, cv.indices = list(), max.threads = 2L, algorithm.config = sgl.standard.config) {
	
	if(length(cv.indices) == 0) {
		
		# Check fold
		if(fold < 2) {
			stop("fold must be equal to or larger than 2")
		}
		
		if(fold > length(data$G)) {
			stop("fold must be equal to or less than the number of samples")
		}
		
		if(fold > max(table(data$G))) {
			# use random sample indices
			use.cv.indices <- TRUE
                        cv.indices <- split(sample(1:(length(data$G))), 1:fold)

                        # TODO need to ensure that each split has one sample from each class

                } else {
                        # compute cv indices
                        cv.indices <- lapply(unique(data$G), function(x) .divide_indices(which(data$G == x), fold))
                        cv.indices <- lapply(cv.indices, function(x) sample(x))
                        cv.indices <- lapply(1:fold, function(i) unlist(lapply(cv.indices, function(x) x[[i]])))
                }

	} else {
          # use user suplied cv splitting
        }

        samples <- 1:max(unlist(cv.indices)) #TODO we should get the number of samples form somewhere else and chek consistency of cv.indices and samples
	training <- lapply(cv.indices, function(x) samples[-x])
	test <- cv.indices
	
	res <- sgl_subsampling(module_name, PACKAGE, data, parameterGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, max.threads, algorithm.config)
	
        # Zero correct and add cv.indices
        res$cv.indices <- cv.indices

        # Reorganize response output
        response_names <- names(res$response[[1]])
        res$response <- lapply(1:length(res$response[[1]]), function(i) lapply(1:length(lambda), function(j) lapply(res$response, function(x) x[[i]][[j]])))

        # Set names
        res$response <- lapply(res$response, function(x) {names(x) <- lambda; x})
        names(res$response) <- response_names

        # Collapse subsample list to matrix
        res$response <- lapply(res$response, function(response_list) .reorder_response_list(response_list, res$cv.indices, data$group.names, data$sample.names))

        # Set class and return
        class(res) <- "sgl"
	return(res)
}

.divide_indices <- function(indices, fold) {

    if(fold > length(indices)) {
        stop("fold large than length of indices vector")
    }

    tmp <- (1:fold)*round(length(indices)/fold)
    tmp[length(tmp)] <- length(indices)
    tmp <- c(0,tmp)

    return(lapply(2:length(tmp), function(i) indices[(tmp[i-1]+1):tmp[i]]))
}

.reorder_response <- function(responses, cv, group.names, sample.names) {

    r <- matrix(nrow = nrow(responses[[1]]), ncol = length(unlist(cv)))
    rownames(r) <- group.names
    colnames(r) <- sample.names

    for(i in 1:length(cv)) {
        r[,cv[[i]]] <- responses[[i]]
    }

    return(r)
}

.reorder_response_list <- function(response_list, cv, group.names, sample.names) lapply(response_list, function(r) .reorder_response(r, cv, group.names, sample.names))
