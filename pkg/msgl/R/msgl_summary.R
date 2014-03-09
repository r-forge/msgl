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

Err <- function(x, type = "count") { #type = count, rate
	
}

features <- function(x) {

  if(class(x) == "msgl") {
      return(sapply(x$beta, function(beta) sum(colSums(beta != 0) != 0)))
  }
  
  stop("Unknown object class")

}

parameters <- function(x) {

  if(class(x) == "msgl") {
      return(sapply(x$beta, function(beta) sum(beta != 0)))
  }
  
  stop("Unknown object class")
}

models <- function(x, index = 1:nmod(x)) {
  return(x$beta[index])
}

coef <- function(x, index = 1:nmod(x)) {

  featurenames = if(is.null(colnames(x$beta[[1]]))) 1:ncol(x$beta[[1]]) else colnames(x$beta[[1]])

  return(lapply(x$beta[index], function(beta) beta[,colSums(beta != 0) != 0]))
}


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

