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

.poissonmodelmoments <- setClass("poissonmodelmoments",
                                 contains = c("dmodelmoments"),
                                 validity = function(object) {
                                     ## else: OK
                                     TRUE
                                 }
)

setMethod("initialize", "poissonmodelmoments",
          function(.Object, ..., model) 
          {
              .Object <- callNextMethod(.Object, ..., model = model)
              generateMoments(.Object)
          }
)

setMethod("generateMoments", "poissonmodelmoments",
          function(object) 
          {
              .generateMomentsPoisson(object)
          } 
)

setMethod("show", "poissonmodelmoments",
          function(object) {
              cat("Object 'modelmoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         :", 
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     factorial   :", 
                  paste(dim(object@factorial), collapse = "x"),
                  "\n")
              cat("     over        :", object@over, "\n")
              cat("     zero        :", object@zero, "\n")
              cat("     model       : Object of class",
                  class(object@model), "\n")
          }
)

## No Setters as users are not intended to manipulate 
## this object ##

### Private functions
### These functions are not exported 
".generateMomentsPoisson" <- function(object) 
{
    K <- object@model@K
    lambda <- object@model@par$lambda
    weight <- object@model@weight
    object@mean <- sum(weight * lambda)
    object@var <- array(sum(weight * lambda * (lambda + 1)) 
                        - object@mean^2, dim = c(1, 1))
    if (K > 1) {
        object@over <- object@var[1] - object@mean 
    } else {
        object@over <- 0
    }
    factm <- array(NA, dim = c(4, 1))
    for (i in seq(1, 4)) {
        factm[i] <- sum(weight * lambda^i)
    }
    dimnames(factm) <- list(c("1st", "2nd", "3rd", "4th"), "")
    object@factorial <- factm
    object@zero <- sum(weight * exp(-lambda))
    return(object)
}

