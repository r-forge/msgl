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

features <- function(object) {

  if(class(object) == "msgl") {
      return(sapply(object$beta, function(beta) sum(colSums(beta != 0) != 0)))
  }
  
  stop("Unknown object class")

}

parameters <- function(object) {

  if(class(object) == "msgl") {
      return(sapply(object$beta, function(beta) sum(beta != 0)))
  }
  
  stop("Unknown object class")
}

nmod <- function(object) {
  return(length(object$beta))
}

models <- function(object, index = 1:nmod(object)) {
  return(object$beta[index])
}

coef <- function(object, index = 1:nmod(object)) {

  featurenames = if(is.null(colnames(object$beta[[1]]))) 1:ncol(object$beta[[1]]) else colnames(object$beta[[1]])

  return(lapply(object$beta[index], function(beta) beta[,colSums(beta != 0) != 0]))
}


print.msgl <- function(object) {

	#FIXME check that  object$msgl_version exsists

  	message(paste("Multinomial logistic regression models (estimated by msgl version ", object$msgl_version, ")", sep=""))
	message()
	message(paste("  This object contains ", nmod(object$beta), " sparse models with ", nrow(object$beta[[1]]), " classes", sep=""))
	message(paste("   between ", min(features(object)), " - ", max(features(object)), " nonzero features and ", min(parameters(object)), " - ", max(parameters(object)), " nonzero parameters.", sep =""))
	message()
	classnames <- paste(rownames(object$beta[[1]])[1:3], collapse=", ") #FIXME what if less than 3 classes + handle long class names
	message(paste("  The first three classes are: ", classnames, sep=""))

}

