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
#' #TODO colSums(res$classes != classes)
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

        #res <- simplify(res)
        #res$link <- lapply(res$link, function(x) t(x))
        #res$response <- lapply(res$response, function(x) t(x))
        #res$classes <- apply(res$classes, MARGIN = 2, function(y) rownames(object$beta[[1]])[y]) #FIXME use classes.names in fit object
        #dimnames(res$classes) <- list(Samples = rownames(x), Lambda.index = 1:length(object$lambda))

        class(res) <- "msgl"
        return(res)
}
