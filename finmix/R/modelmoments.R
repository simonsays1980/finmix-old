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
         representation(mean    = "vector",
                        var     = "array",
                        model   = "model"
                        ),
         validity = function(object) {
             ## else: OK
		     TRUE
         },
         prototype(
                   mean     = vector(),
                   var      = array(),
                   model    = model()
                   )
)

"modelmoments" <- function(model) {
    dist <- model@dist
    if (dist == "normult") {
        object <- .normultmodelmoments(model = model)
    } else if (dist == "studmult") {
        object <- .studmultmodelmoments(model = model)
    } else if (dist == "student") {
        object <- .studentmodelmoments(model = model)
    } else if (dist == "normal") {
        object <- .normalmodelmoments(model = model)
    } else if (dist == "exponential") {
        object <- .exponentialmodelmoments(model = model)
    }
    return(object)
}

## Getters ##
setGeneric("getMean", function(object) standardGeneric("getMean"))
setMethod("getMean", "modelmoments", 
          function(object) {
              return(object@mean)
          }		
)
setGeneric("getVar", function(object) standardGeneric("getVar"))
setMethod("getVar", "modelmoments", 
          function(object) {
              return(object@var)
          }
)

setGeneric("getModel", function(object) standardGeneric("getModel"))
setMethod("getModel", "modelmoments",
          function(object) {
              return(object@model)
          }
)
## Setters are not provided as users are not intended to manipulate ##
## this object ##
