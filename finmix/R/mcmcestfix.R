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

setClass("mcmcestfix",
         representation(dist        = "character",
                        K           = "integer",
                        indicmod    = "character",
                        map         = "list",
                        bml         = "list",
                        ieavg       = "list"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("show", "mcmcestfix", 
          function(object) {
              cat("Object 'mcmcest\n")
              cat("     dist        :", object@dist, "\n")
              cat("     K           :", object@K, "\n")
              cat("     indicmod    :", object@indicmod, 
                  "\n")
              cat("     map         : List of", 
                  length(object@map), "\n")
              cat("     bml         : List of",
                  length(object@bml), "\n")
              cat("     ieavg       : List of", 
                  length(object@ieavg), "\n")
          }
)

## Getters ##
## Generic set in 'model.R' ##
setMethod("getDist", "mcmcestfix", 
          function(object) {
              return(object@dist)
          }
)
## Generic set in 'model.R' ##
setMethod("getK", "mcmcestfix", 
          function(object) {
              return(object@K)
          }
)
## Generic set in 'model.R'
setMethod("getIndicmod", "mcmcestfix", 
          function(object) {
              return(object@indicmod)
          }
)

setGeneric("getMap", function(object) standardGeneric("getMap"))
setMethod("getMap", "mcmcestfix", 
           function(object) {
               return(object@map)
           }
)

setGeneric("getBml", function(object) standardGeneric("getBml"))
setMethod("getBml", "mcmcestfix",
          function(object) {
              return(object@bml)
          }
)

setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))
setMethod("getIeavg", "mcmcestfix", 
          function(object) {
              return(object@ieavg)
          }
)

## No setters as users are not intended to manipulate
## this object

