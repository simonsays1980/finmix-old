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

## 'datamoments' is a virtual class from which the corresponding      ##
## datamoments for 'continuous' and 'discrete' inherit	 		## 
setClass("datamoments",
         representation(mean            = "numeric",
                        var             = "matrix",
		                data            = "data",
                		"VIRTUAL")
)

## Getters ##
## Generic set in 'groupmoments.R' ##
setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Setters ##
## No setters as users should not manipulate a 'datamoments' object ##

## mutual constructor for all type of datamoments ##
"datamoments" <- function(value = data()) {
        if (all(is.na(value@y))) {
            stop("'data' object has no data. Slot 'y' is empty.")
        }
		if (value@type == "continuous") {
			.Object <- .cdatamoments(value = value)
        } else { 
			.Object <- .ddatamoments(value = value)
        }
		return(.Object)
}


