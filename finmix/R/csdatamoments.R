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

setClass("csdatamoments",
         representation(B       = "vector",
		                W       = "vector",
		                T       = "vector",
                        R       = "numeric",
		                Rtr     = "numeric",
		                Rdet    = "numeric"),
         contains = c("sdatamoments"),
	     validity = function(object) {
             ## else: ok
		     TRUE
         }
)

setMethod("initialize", "csdatamoments",
          function(.Object, value) {
              object <- as(.Object, "sdatamoments")
              object <- callNextMethod(object, value)
              as(.Object, "sdatamoments") <- object
              .Object <- generateMoments(.Object)
              return(.Object)
          }
)

setMethod("generateMoments", "csdatamoments",
          function(object) {
	          ## enforce column.wise ordering ##
              if (!object@data@bycolumn) {
		          datam     <- t(object@data@y)
		          classm    <- t(object@data@S)
              } else {
        		  datam     <- object@data@y
		          classm    <- object@data@S
	          }
   	          ## Calculate the between-group variance ##
	          ## 'B' is an r x r matrix ##
              gmeans <- object@gmoments@mean
              nkm <- object@gmoments@NK
              ## Calculate the total heterogeneity ##
              ## 'T' is an r x r array ##
              object@T <- var(datam, na.rm = TRUE) * nrow(datam)
  		      ## Calculate the within-group heterogeneity ##
			  ## 'W' is an r x r array ##
              wkm <- object@gmoments@WK
              object@W <- apply(wkm, c(1, 2), sum, na.rm = TRUE)
              ## Calculate between-group heterogeneity ##
			  ## 'B' is an r x r array ##
			  object@B <- object@T - object@W
			  ## Calculate coefficient of determination ##
			  ## 'Rtr' is an 1 x 1 numeric ##
			  ## 'Rdet' is an 1 x 1 numeric ##
			  if (object@data@r > 1) {
                  r <- NA
                  object@R <- as.numeric(r)
                  object@Rtr <- 1 - sum(diag(object@W), na.rm = TRUE) / 
                                    sum(diag(object@T), na.rm = TRUE) 
			      object@Rdet <- 1 - det(object@W)/det(object@T)
              } else {
                  rtr <- NA
                  rdet <- NA
                  object@Rtr <- as.numeric(rtr)
                  object@Rdet <- as.numeric(rdet)
                  object@R <- 1 - object@W[1]/object@T[1] 
              }
			  return(object)
          }
)

setMethod("show", "csdatamoments", 
          function(object) {
              cat("Object 'sdatamoments'\n")
              cat("     class       :", class(object), "\n")
              cat("     B           : Vector of", 
                  length(object@B), "\n")
              cat("     W           : Vector of",
                  length(object@W), "\n")
              cat("     T           : Vector of",
                  length(object@T), "\n")
              if (object@data@r > 1) {
                  cat("     Rdet        :", object@Rdet, "\n")
                  cat("     Rtr         :", object@Rtr, "\n")
              }
              cat("     gmoments    : Object of class", 
                  class(object@gmoments), "\n")
              cat("     data        : Object of class", 
                  class(object@data), "\n")
          }
)

## Getters ##
setMethod("getGmoments", "csdatamoments", function(object) {
						return(object@gmoments)
					}
)
setMethod("getWK", "csdatamoments", function(object) {
						return(object@WK)
					}
)
setMethod("getVar", "csdatamoments", function(object) {
				return(object@var)
			}
)
setGeneric("getB", function(object) standardGeneric("getB"))
setMethod("getB", "csdatamoments", function(object) {
						return(object@B)
					}
)
setGeneric("getW", function(object) standardGeneric("getW"))
setMethod("getW", "csdatamoments", function(object) {
						return(object@W)
					}
)
## Already set as generic in 'model.R' ## 
setMethod("getT", "csdatamoments", function(object) {
						return(object@T)
					}
)
## Generic set in 'model.R' ##
setMethod("getR", "csdatamoments",
          function(object) {
              return(object@R)
          }
)
setGeneric("getRtr", function(object) standardGeneric("getRtr"))
setMethod("getRtr", "csdatamoments", function(object) {
						return(object@Rtr)
					}
)
setGeneric("getRdet", function(object) standardGeneric("getRdet"))
setMethod("getRdet", "csdatamoments", function(object) {
						return(object@Rdet)
					}
)

## Generic set in 'groupmoments.R' ##
setMethod("getData", "csdatamoments", function(object) {
						return(object@data)				
					}
)

## Setters ##
## No setters, as it users are adviced not to manipulate moment objects ##
