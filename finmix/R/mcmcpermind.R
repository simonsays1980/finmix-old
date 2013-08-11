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

setClass("mcmcpermind", 
         representation(weightperm  = "array",
                        entropyperm = "array",
                        STperm      = "array",
                        Sperm       = "array",
                        NKperm      = "array"
                        ),
         contains = c("mcmcpermfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

## Getters ##
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))
setMethod("getWeightperm", "mcmcpermind", function(object) {
          return(object@weightperm)
})

setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))
setMethod("getEntropyperm", "mcmcpermind", 
          function(object) {
              return(object@entropyperm)
          }
)

setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))
setMethod("getSTperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

setGeneric("getSperm", function(object) standardGeneric("getSperm"))
setMethod("getSperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))
setMethod("getNKperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
