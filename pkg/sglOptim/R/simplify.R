#
#     Description of this R script:
#     Routines for handling simplifying responses
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

#TODO dim names

#TODO generic s3 function
#sgl_simplify <- function(x, ...) UseMethod("sgl_simplify")


sgl_simplify_list <- function(x, ...) {

    require(plyr)

    return(laply(x, .fun = function(y) y))
}

simplify <- function(x) {

    require(plyr)

    if(class(x[[1,1]]) == "list") {
        response.names <- names(x[[1,1]])

        res <- llply(response.names, function(name)  simplify(llply(x, .fun = function(y) y[[name]])))
        names(res) <- response.names

        return(res)
    }

    if(is.vector(class(x[[1,1]])) && length(x[[1,1]]) > 1) {
        return(alply(x, .margins = 2, .fun = sgl_simplify_list))
    }

    if(is.vector(class(x[[1,1]]))) {
        return(aaply(x, .margins=c(1,2), .fun = function(y) y))
    }

    stop("Cannot simplify object")
}

