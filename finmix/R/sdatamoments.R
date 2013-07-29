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
	representation(
		groupmoments = "groupmoments",
		WK = "array",
		var = "array",
		data = "data"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

## mutual constructor for all type of datamoments ##
"sdatamoments" <- function(data) {
			
			if(data@type == "discrete") {
				groupmoments <- groupmoments(data)
				## Calculate moments ##
				## work only with data ordered by column ##
				if(!data@bycolumn) {
					datam <- t(data@y)
					classm <- t(data@S)
				}
				else {
					datam <- data@y
					classm <- data@S				
				}
				has.colnames <- !is.null(colnames(datam))
				## Calculate within-group variability ##
				## 'WK' is an 1 x K array for r == 1 ##
				## 'WK' is an r x r x K array for r > 1 ##
				## Calculate the variance ##
				## 'var' is an 1 x K array for r == 1 ##
				## 'var' is an r x r x K array for r > 1 ##
				level.set <- as.numeric(levels(factor(classm)))
				K <- length(level.set)
				kmeans <- groupmoments@group.mean
				r <- ncol(datam)
				nkm <- groupmoments@NK
				if(r == 1) {
					wkm <- array(0, dim = c(1, K))
					varm <- array(0, dim = c(1, K))
					for(i in 1:K) {
						group <- matrix(datam[which(classm == level.set[i]),])
						wkm[i] <- sum((group - kmeans[i])^2, na.rm = TRUE)
						varm[i] <- wkm[i]/kmeans[i]
					}
					dimnames(wkm) <- list(NULL, paste("k=", 1:K, sep=""))
					dimnames(varm) <- list(NULL, paste("k=", 1:K, sep=""))
					if(has.colnames) {
						rownames(wkm) <- colnames(datam)
						rownames(varm) <- colnames(datam)
					}
				}
				else { ## r > 1
					wkm <- array(0, c(r, r, K))
					varm <- array(0, c(r, r, K))
					for(i in 1:K) {
						kmeanm <- t(as.matrix(kmeans[,i])[, rep(1, nkm[i])])
						group <- as.matrix(datam[which(classm == level.set[i]),]) - kmeanm
						out.prod <- matrix(0, ncol  = r, nrow = r)
						for(j in 1:r) {
							out.prod <- out.prod + group[j,] %o% group[j,] 	
						}
						wkm[,,i] <- out.prod
						varm[,,i] <- out.prod/nkm[i] 
					}
				}
				.Object <- new("sdatamoments", groupmoments = groupmoments, WK = wkm, 
						var = varm, data = data)
			}
			else { ## type = "continuous"
				.Object <- .csdatamoments(data)
			}
			return(.Object)
}

setMethod("show", "sdatamoments", function(object) {
                                                name.data <- getName(object@data)
                                                oname <- ifelse(length(name.data) == 0, "", name.data)
                                                cat("SDataMoments object '", oname, "'\n")
                                                gmoments <- object@groupmoments
                                                cat("   Type            :", class(object), "\n")
                                                cat("   Group Moments   :", class(object@groupmoments), "\n")
                                                cat("   WK              :", paste(dim(object@WK), collapse="x"), "\n")
                                                cat("   Var (within)    :", paste(dim(object@var), collapse="x"), "\n")                                         }
)

## Getters ##
setGeneric("getGroupMoments", function(.Object) standardGeneric("getGroupMoments"))
setMethod("getGroupMoments", "sdatamoments", function(.Object) {
							return(.Object@groupmoments)
						}
)
setGeneric("getWK", function(.Object) standardGeneric("getWK"))
setMethod("getWK", "sdatamoments", function(.Object) {
						return(.Object@WK)
					}
)
## Generic set in 'modelmoments' class ##
setMethod("getVar", "sdatamoments", function(.Object) {
						return(.Object@var)
					}
)

## Setters ##
## No Setters, as it is adviced for users not to manipulate moment objects ##
