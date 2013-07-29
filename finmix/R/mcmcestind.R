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

setClass("mcmcestind",
         representation(eavg = "list"),
         contains = c("mcmcestfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setClassUnion("mcmcest", c("mcmcestfix",
                           "mcmcestind")
)

setMethod("show", "mcmcestind",
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
              cat("     eavg        : List of",
                  length(object@eavg), "\n")
          }
)

## Getters ##
setGeneric("getEavg", function(object) standardGeneric("getEavg"))
setMethod("getEavg", "mcmcestind", 
          function(object) {
              return(object@eavg)
          }
)

## No setters as users are not intended to manipulate 
## this object.
