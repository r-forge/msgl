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
#' @param collapse if \code{TRUE} the results for each subsample will be collapse into one result (this is useful if the subsamples are not overlapping)
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
#' test <- replicate(2, sample(1:length(classes))[1:20], simplify = FALSE)
#' train <- lapply(test, function(s) (1:length(classes))[-s])
#'
#' fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda,
#'  training = train, test = test)
#'
#' # Missclassification count of second subsample
#' # subsample 1
#' colSums(fit.sub$classes[[1]] != classes[test[[1]]])
#' # subsample 2
#' colSums(fit.sub$classes[[2]] != classes[test[[2]]])
#' @author Martin Vincent
#' @export
#' @useDynLib msgl .registration=TRUE
msgl.subsampling <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE,
                lambda, training, test, sparse.data = is(x, "sparseMatrix"), collapse = FALSE, max.threads = 2L, algorithm.config = sgl.standard.config) {
	
        # Default values
        if(is.null(grouping)) {
                covariateGrouping <- factor(1:ncol(x))
        } else {
                # ensure factor
                covariateGrouping <- factor(grouping)
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
        data <- create.sgldata(x, y = NULL, sampleWeights, classes, sparseX = sparse.data)

        # call sglOptim function
        if(data$sparseX) {

                res <- sgl_subsampling("msgl_sparse", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, collapse, max.threads, algorithm.config)
        } else {

                res <- sgl_subsampling("msgl_dense", "msgl", data, covariateGrouping, groupWeights, parameterWeights, alpha, lambda, training, test, collapse, max.threads, algorithm.config)
        }

		### Responses
		res$classes <- res$responses$classes
		res$response <- res$responses$response
		res$link <- res$responses$link
		res$responses <- NULL
		
		# Set class names
		if(!is.null(data$group.names)) {
			for(i in 1:length(training)) {
				res$classes[[i]] <- apply(X = res$classes[[i]], MARGIN = c(1,2), FUN = function(x) data$group.names[x])
				res$link[[i]] <- lapply(X = res$link[[i]], FUN = function(m) {rownames(m) <- data$group.names; m})
				res$response[[i]] <- lapply(X = res$response[[i]], FUN = function(m) {rownames(m) <- data$group.names; m})
			}
		}
		
        res$msgl_version <- packageVersion("msgl")

        class(res) <- "msgl"
        return(res)
}
