## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

## 'datamoments' is a virtual classes from which the corresponding      ##
## datamoments for 'continuous' and 'discrete' inherit	 		## 
setClass("sdatamoments",
         representation(gmoments    = "groupmoments",
                        data        = "data"),
         validity = function(object) {
             ## else: OK
		     TRUE
         }
)

## mutual constructor for both types of sdatamoments ##
"sdatamoments" <- function(value = data()) 
{
    if (all(is.na(value@y))) {
        stop("'data' object has no data. Slot 'y' is empty.")
    } else {
        if (all(is.na(value@S))) {
            stop("'data' object has no allocations. Slot 'S' is empty.")
        }
    }
    if (value@type == "discrete") {
        object <- new("sdatamoments", value = value)
    } else {
        object <- .csdatamoments(value = value)
    }
    return(object)
}

setMethod("initialize", "sdatamoments",
          function(.Object, ..., value = data()) 
          {
              .Object@data <- value
              .Object@gmoments <- new("groupmoments", value = value)
              return(.Object)
          }
)

setMethod("show", "sdatamoments", 
          function(object) 
          {
              cat("Object 'sdatamoments'\n")
              cat("     gmoments    : Object of class",
                  class(object@gmoments), "\n")
              cat("     data        : Object of class",
                  class(object@data), "\n")
          }
)

## Getters ##
setMethod("getGmoments", "sdatamoments", function(object) {
							return(object@gmoments)
						}
)

setMethod("getData", "sdatamoments", 
          function(object) {
              return(object@data)
          }
)

## Setters ##
## No Setters, as it is adviced for users not to manipulate moment objects ##
