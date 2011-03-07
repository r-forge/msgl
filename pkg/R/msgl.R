# TODO: Add comment
# 
# Author: martin
###############################################################################


dyn.load("/home/martin/SGL/SGLRInterface/Release/libSGLRInterface.so")

sparseGroupLasso.experimental <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, refit = FALSE, delta = 0.001) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	betaList <- .Call("r_sgl_experimental", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, refit)
	
	for(l in 1:numberOfModels) {
		
		if(standardize) {
			
			beta.org <- betaList[[l]] %*% diag(c(1,1/x.scale))
			beta.org[,1] <- beta.org[,1] - rowSums(betaList[[l]][,-1] %*% diag(x.center/x.scale))
			
			betaList[[l]] <- beta.org
			
		}
		
		colnames(betaList[[l]]) <- c("(Intercept)", colnames(x));
		rownames(betaList[[l]]) <- levels(classes);
		
	}
	
	return(betaList)
}

sparseGroupLasso.experimental.cv <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, delta = 0.001, refit = FALSE, fold = 10L, numberOfThreads = 2L) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	predictedClasses <- .Call("r_sgl_experimental_cv", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, refit, fold, numberOfThreads)
	
	return(predictedClasses)
}

sparseGroupLasso <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, refit = FALSE, delta = 0.001) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	betaList <- .Call("r_sgl_simple", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, refit)
	
	for(l in 1:numberOfModels) {
		
		if(standardize) {
			
			beta.org <- betaList[[l]] %*% diag(c(1,1/x.scale))
			beta.org[,1] <- beta.org[,1] - rowSums(betaList[[l]][,-1] %*% diag(x.center/x.scale))
			
			betaList[[l]] <- beta.org
			
		}
		
		colnames(betaList[[l]]) <- c("(Intercept)", colnames(x));
		rownames(betaList[[l]]) <- levels(classes);
		
	}
	
	return(betaList)
}

sparseGroupLasso.cv <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, delta = 0.001, refit = FALSE, fold = 10L, numberOfThreads = 2L) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	predictedClasses <- .Call("r_sgl_simple_cv", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, refit, fold, numberOfThreads)
	
	return(predictedClasses)
}

sparseGroupLasso.stability.paths <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, delta = 0.001, numberOfSubsets = 100L, numberOfThreads = 2L) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	paths <- .Call("r_sgl_simple_stability_paths", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, numberOfSubsets, numberOfThreads)
	
	for(l in 1:numberOfModels) {
		colnames(paths[[l]]) <- c("(Intercept)", colnames(x));
		rownames(paths[[l]]) <- levels(classes);
	}
	
	return(paths)
}

sparseGroupLasso.stability.selection <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, delta = 0.001, cutoff = 0.7, numberOfSubsets = 100L, numberOfThreads = 2L) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	betas <- .Call("r_sgl_simple_stability_selection", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, cutoff, numberOfSubsets, numberOfThreads)
	
	for(l in 1:numberOfModels) {
		
		if(standardize) {
			
			beta.org <- betas[[l]] %*% diag(c(1,1/x.scale))
			beta.org[,1] <- beta.org[,1] - rowSums(betas[[l]][,-1] %*% diag(x.center/x.scale))
			
			betas[[l]] <- beta.org
			
		}
		
		colnames(betas[[l]]) <- c("(Intercept)", colnames(x));
		rownames(betas[[l]]) <- levels(classes);
		
	}
	
	return(betas)
}

sparseGroupLasso.stability.selection.cv <- function(x, classes, lambdaMin = 0.02, numberOfModels = 50L, alpha = 0.5, featureWeights = rep(1, ncol(x)), classWeights = rep(1, length(levels(classes))), standardize = TRUE, delta = 0.001, fold = 10L, cutoff = 0.7, numberOfSubsets = 100L, numberOfThreads = 2L) {
	
	if(!is.factor(classes)) {
		stop("Parameter classes must be a factor - ")
	}
	
	if(standardize) {
		x <- scale(x, TRUE, TRUE)
		x.scale <- attr(x,"scaled:scale")
		x.center <- attr(x,"scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		stop("NA values in x")
	}
	
	predictedClasses <- .Call("r_sgl_simple_stability_selection_cv", x, as.integer(classes)-1L, numberOfModels, lambdaMin, alpha, featureWeights, classWeights, delta, fold, cutoff, numberOfSubsets, numberOfThreads)
	
	return(predictedClasses)
}

predict.sgl.test <- function(beta.list, x) {
	
	numberOfClasses <- dim(beta.list[[1]])[1]
		
	prop <- lapply(1:length(beta.list), function(i) .Call("r_sgl_predict_classes", x, beta.list[[i]]))
	return(prop)
}

#' Computes the predicted probabilities/classes
#' @param beta.list 
#' @param x 
#' @param type 
#' @return a list
#' @author MVI
#' @export
predict.sgl <- function(beta.list, x, type = "class") {
	
	#Intercept
	x <- cbind(rep(1, nrow(x)), x)
	colnames(x)[1] <- "(Intercept)"
	
	numberOfClasses <- dim(beta.list[[1]])[1]
	classNames <- rownames(beta.list[[1]])
	
	if(type == "response") {
		prop <- lapply(1:length(beta.list), function(l) sapply(1:numberOfClasses, function(i) prop(x %*% t(beta.list[[l]]), i)))
		return(prop)
	}
	
	if(type == "link") {
		link <- lapply(1:length(beta.list), function(l) sapply(1:numberOfClasses, function(i) (x %*% t(beta.list[[l]]))[,i] ))
		return(link)
	}
	
	if(type == "class") {
		prop <- lapply(1:length(beta.list), function(l) sapply(1:numberOfClasses, function(i) prop(x %*% t(beta.list[[l]]), i)))
		lapply(1:length(prop), function(l) sapply(1:nrow(prop[[l]]), function(i) classNames[.argmax(prop[[l]][i,])[1]]))
	}
}

.argmax <- function(x) which(max(x) == x)

#' Calculate the partition function
#' @param eta a matrix. 
#' @return a vector of length equal to the rowdim of eta. The i'th element of the vector is the value of the partition function on the i'th row of eta.
#' @author MVI
#' @export
z <- function(eta) rowSums(exp(eta))

#' Calculate the estimated probability of the given class
#' @param eta a matrix
#' @param class
#' @return a vector of probabilities 
#' @author MVI
#' @export
prop <- function(eta, class) exp(eta[,class])/z(eta)
