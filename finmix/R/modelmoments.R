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

setClass("modelmoments", 
	representation(
	mean = "matrix",
	var = "matrix"
	),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

"modelmoments" <- function(model, J. = 4) {
		dist <- model@dist
		if(dist == "normal" || dist == "normult" ) {
			
			nsmodelmoments <- nsmodelmoments(model, J.)
			return(nsmodelmoments)
		}
		else if(dist == "student" || dist == "studmult") {
			return("[Warning] Function 'modelmoments' not implemented for Student mixtures.")
		}
		else if(dist == "poisson") {
			
			dmodelmoments <- dmodelmoments(model, J.)
			return(dmodelmoments)
		}
		else if(dist == "binomial") {
			return("[Warning] Function 'modelmoments' not implemented for binomial mixtures.")
		}
		else if(dist == "exponential") {
			return("[Warning] Function 'modelmoments' not implemented for exponential mixtures")
		}
		
}

## Getters ##
setGeneric("getMean", function(.Object) standardGeneric("getMean"))
setMethod("getMean", "modelmoments", function(.Object) {
						return(.Object@mean)
					}		
)
setGeneric("getVar", function(.Object) standardGeneric("getVar"))
setMethod("getVar", "modelmoments", function(.Object) {
						return(.Object@var)
					}
)

## Setters are not provided as users should not manipulate a 'modelmoments' object ##
