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

simplify_vector <- function(responses) {

    require(plyr)

    return(laply(responses, .fun = function(x) x))
}

simplify <- function(res, simplify.fun) {

    require(plyr)

    return(alply(res, .margins = 2, .fun = simplify.fun))
}
