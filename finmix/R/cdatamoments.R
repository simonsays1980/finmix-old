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

setClass("cdatamoments",
	representation(
		higher = "matrix",
		skewness = "matrix",
		kurtosis = "matrix",
		corr = "matrix"		
		),
	contains = c("datamoments"),
	validity = function(object) {
		mom.moments <- object@higher
		mom.data <- object@data
		mom.r <- mom.data@r
		if(mom.r != ncol(object@mean)) 
			return("[Error] Data dimension and dimension of the mean do not match.")
		if(mom.r != ncol(object@var)) 
			return("[Error] Data dimension and dimension of the variance do not match.")
		if(mom.r != ncol(mom.moments) || nrow(mom.moments) != 4) 
			return("[Error] Data dimension and dimension of the L=4 moments do not match.")
		if(mom.r != ncol(object@skewness)) 
			return("[Error] Data dimension and dimension of the skewness do not match.")
		if(mom.r != ncol(object@kurtosis)) 
			return("[Error] Data dimension and dimension of the kurtosis do not match.")
		if(!is.na(object@corr) && mom.r == 1)
			return("[Error] Data has dimension 1. No correlation matrix available.")
		if(!is.na(object@corr) && (mom.r != ncol(object@corr) || mom.r != nrow(object@corr)))
			return("[Error] Data dimension and dimension of correlation matrix do not match.")
		if(!is.na(object@kurtosis) && any(object@kurtosis < 0))
			return("[Error] Kurtosis is negative.")
		if(!is.na(object@var) && any(object@var < 0))
			return("[Error] Variance is negative")
		if(!is.na(mom.moments) && any(mom.moments[2,] < 0))
			return("[Error] Second higher moment is negative.")
		if(!is.na(mom.moments) && any(mom.moments[4,] < 0))
			return("[Error ]Fourth higher moment is negative.")
		## else: ok
		TRUE
	}
)

setMethod("initialize", "cdatamoments", function(.Object, ..., data) {
						.Object@data <- data
						cmoments(.Object) <- .Object@data
						has.S <- !all(is.na(data@S))
						if(has.S) {
							.Object@sdatamoments <- sdatamoments(data)
						}
						callNextMethod(.Object, ...)
						return(.Object)
					}
)

"cmoments<-" <-  function(.Object, value) {
			              
				## Compute means ##
				## work only with data ordered by column ##
				if(!value@bycolumn) {
					datam <- t(value@y)
				}
				else {
					datam <- value@y
				}
				has.colnames <- !is.null(colnames(datam))
				## means is a matrix r x 1 ##
				means <- matrix(0, nrow = value@r, ncol = 1)
				for(i in 1:value@r) {
					means[i] <- mean(datam[,i])
				}
				if(has.colnames) {
					rownames(means) <- colnames(datam)
				}
				.Object@mean <- means
				
				## Compute variances or variance-covariance matrix ##
				## variance is a r x r matrix ##
				
				.Object@var <- t(as.matrix(diag(var(datam))))
				
				## Compute higher moments ##
				## higher.moments is a r x L matrix (L = 4) ##
				d = datam - mapply(rep, means, nrow(datam)) 
				momentsm <- matrix(0, nrow = value@r, ncol = 4)
				for(i in 1:getR(value)) {
					momentsm[i, 2] <- mean(d[,i]^2)
					momentsm[i, 3] <- mean(d[,i]^3)
					momentsm[i, 4] <- mean(d[,i]^4) 
				}
				if(has.colnames) {
					rownames(momentsm) <- colnames(datam)
				}
				.Object@higher <- momentsm

				## Compute skewness and kurtosis ##
				## skewness and kurtosis are r x 1 matrices ## 
				skewm <- matrix(0, nrow = value@r, ncol = 1)
				kurtm <- matrix(0, nrow = value@r, ncol = 1)
				for(i in 1:value@r) {
					skewm[i] <- momentsm[i, 3]/momentsm[i, 2]^1.5
					kurtm[i] <- momentsm[i, 4]/momentsm[i, 2]^2
				}
				if(has.colnames) {
					rownames(skewm) <- colnames(datam)
					rownames(kurtm) <- colnames(datam) 	
				}
				.Object@skewness <- skewm
				.Object@kurtosis <- kurtm

				## Compute corr matrix in case of r > 1 ##
				## corr is a r x r matrix ##
				if(value@r > 1) {
					.Object@corr <- cor(value@y)			
				}
				else {
					.Object@corr <- matrix()
				}

				return(.Object)
}

setMethod("show", "cdatamoments", 
          function(object) {
		      name.data <- getName(getData(object))
			  oname <- ifelse(length(name.data) == 0, "datamoments", name.data)
              cat("Object '", oname, "'\n")
              cat("     class       :", class(object), "\n")
              cat("     mean        :", 
                  paste(dim(object@mean), collapse = "x"), "\n")
              cat("     var         :", 
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     higher      :", 
                  paste(dim(object@higher), collapse = "x"), "\n")
              cat("     skewness    :", 
                  paste(dim(object@skewness), collapse = "x"), "\n")
              cat("     kurtosis    :", 
                  paste(dim(object@kurtosis), collapse = "x"), "\n")
              cat("     corr        :", 
                  paste(dim(object@corr), collapse = "x"), "\n")
              if (!all(is.na(object@data@S))) {
                  cat("     sdatamoments    :", class(object@sdatamoments), "\n")
              }
          }
)

## Getters ##
## Generic set in 'modelmoments' class ##
setMethod("getMean", "cdatamoments", function(.Object) {
						return(.Object@mean)
					}
)
## Generic set in 'modelmoments' class ##
setMethod("getVar", "cdatamoments", function(.Object) {
						return(.Object@var)
					}
)
## Generic set in 'datamoments' class ##
setMethod("getData", "cdatamoments", function(.Object) {
				return(.Object@data)
			}
)
## Generic set in 'datamoments' class ##
setMethod("getSDatamoments", "cdatamoments", function(.Object) {
							return(.Object@sdatamoments)
						}
)

## Generic set in 'nsmodelmoments' class ##
setMethod("getHigher", "cdatamoments", function(.Object) {
						return(.Object@higher)
					}
)
## Generic set in 'nsmodelmoments' class ##
setMethod("getSkewness", "cdatamoments", function(.Object) {
						return(.Object@skewness)
					}
)
## Generic set in 'nsmodelmoments' class ##
setMethod("getKurtosis", "cdatamoments", function(.Object) {
						return(.Object@kurtosis)
					}
)
## Generic set in 'nsmodelmoments' class ##
setMethod("getCorr", "cdatamoments", function(.Object) {
						return(.Object@corr)
					}
)
## Setters ##
## No setters as users should not manipulate a 'cdatamoments' object ##
