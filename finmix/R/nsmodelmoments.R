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

setClass("nsmodelmoments",
	representation(
	higher = "matrix",
	skewness = "matrix",
	kurtosis = "matrix",
	R = "numeric",
	Rtr = "numeric",
	Rdet = "numeric",
	corr = "matrix"),
	contains = c("modelmoments"),
	validity = function(object) {
		## else: OK
		TRUE
	}	
)

"nsmodelmoments" <- function(model, J.) {
			dist <- model@dist
			
			if(dist == "normal") {
				
                       		meanm <- matrix(sum(model@weight * model@par$mu))
                        	higher <- mixturemoments.normal(model, J., meanm)
                        	varm <- matrix(higher[2])
					
				skewness <- matrix()
				kurtosis <- matrix()
				if(J. > 2) {
					skewness <- matrix(higher[3]/higher[2]^1.5)
				}
				if(J. > 3) {
					kurtosis <- matrix(higher[4]/higher[2]^2)
				}
				B <- sum(model@weight * (model@par$mu - matrix(meanm, nrow = 1, ncol = model@K))^2)
				W <- sum(model@weight * model@par$sigma)
				R <- 1 - W/varm[1]
				
				nsmodelmoments <- new("nsmodelmoments", mean = meanm, var = varm, higher = higher
							, skewness = skewness, kurtosis = kurtosis, R = R
							, Rtr = numeric(), Rdet = numeric(), corr = matrix()) 
			}
			
			else if(dist == "normult") {
				
				r <- model@r
				meanm <- apply(model@par$mu, 1, "*", model@weight)
				meanm <- matrix(apply(meanm, 1, sum))
				varm <- matrix(0, nrow = r, ncol = r)
				W <- matrix(0, nrow = r, ncol = r)
				B <- matrix(0, nrow = r, ncol = r)
				
				for(k in 1:model@K) {
					varm <- varm + model@par$mu[,k] %*% t(model@par$mu[,k]) 
						+ model@par$sigma[,,k] * model@weight[k] 
					W <- W + model@par$sigma[,,k] * model@weight[k]
					d <- model@par$mu[,k] - meanm
					B <- B + d %*% t(d) * model@weight[k]  
				} 
				varm <- varm - meanm %*% t(meanm)
				
				cd <- diag(1/diag(varm)^.5)
				corr <- cd %*% varm %*% cd

				Rtr <- 1 - sum(diag(W))/sum(diag(varm))
				Rdet <- 1 - det(W)/det(varm)

				higher <- matrix(0, nrow = r, ncol = J.)
				
				if(J. < 3) { ## set skewness and kurtosis to 'NA'
					skewness <- matrix()
					kurtosis <- matrix()

					for(i in 1:model@r) {
						marmodel <- mixturemar(model, i) 
						higher[i,] <- t(mixturemoments.normal(marmodel, J., meanm[i,]))
					}
				}
				if(J. == 3) { ## set kurtosis to 'NA'
					skewness <- matrix(0, nrow = r, ncol = 1)
					kurtosis <- matrix()
					
					for(i in 1:model@r) {
						marmodel <- mixturemar(model, i)
						higher[i,] <- t(mixturemoments.normal(marmodel, J., meanm[i,]))
						skewness[i] <- higher[i,3]/higher[i,2]^1.5
					}
				}
				else {
					skewness <- matrix(0, nrow = r, ncol = 1)	
					kurtosis <- matrix(0, nrow = r, ncol = 1)
				
					for(i in 1:model@r) {
						marmodel <- mixturemar(model, i)
						higher[i,] <- t(mixturemoments.normal(marmodel, J., meanm[i,]))
						skewness[i] <- higher[i,3]/higher[i,2]^1.5
						kurtosis[i] <- higher[i,4]/higher[i,2]^2
					}
				}
									
				nsmodelmoments <- new("nsmodelmoments", mean = meanm, var = varm, 
							higher = higher, skewness = skewness, kurtosis = kurtosis,
							R = numeric(), Rtr = Rtr, Rdet = Rdet, corr = corr)	
				
			} 
}
	
setMethod("show", "nsmodelmoments", function(object) {
						cat("Modelmoments object\n")
						cat("	Type	:", class(object), "\n")
						cat("	Mean	:", paste(dim(object@mean), collapse = "x"), "\n")
						cat("	Var	:", paste(dim(object@var), collapse = "x"), "\n")
						cat("	Higher	:", paste(dim(object@higher), collapse = "x"), "\n")
						if(length(object@R) > 0) {
							cat("	R	: [", object@R, "]\n")
						}
						else {
							cat("	Rtr	: [", object@Rtr, "]\n")
							cat("	Rdet	: [", object@Rdet, "]\n")
							cat("	Corr	:", paste(dim(object@corr), collapse = "x"), "\n")
						}
					}
)
## Getters ##
## Generic set in 'modelmoments' class ##
setMethod("getMean", "nsmodelmoments", function(object) {
						return(object@mean)
					}
)
## Generic set in 'modelmoments' class ##
setMethod("getVar", "nsmodelmoments", function(object) {
						return(object@var)
					}
)
setGeneric("getHigher", function(object) standardGeneric("getHigher"))
setMethod("getHigher", "nsmodelmoments", function(object) {
							return(object@higher)
						}
)
setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))
setMethod("getSkewness", "nsmodelmoments", function(object) {
							return(object@skewness) 
						}
)
setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))
setMethod("getKurtosis", "nsmodelmoments", function(object) {
							return(object@kurtosis)
						}
)
setGeneric("getCorr", function(object) standardGeneric("getCorr"))
setMethod("getCorr", "nsmodelmoments", function(object) {
						return(object@corr)
					}
)
setGeneric("getRdet", function(object) standardGeneric("getRdet"))
setMethod("getRdet", "nsmodelmoments", 
          function(object) {
              return(object@Rdet)
          }
)
setGeneric("getRtr", function(object) standardGeneric("getRtr"))
setMethod("getRtr", "nsmodelmoments", 
          function(object) {
              return(object@Rtr)
          }
)
## Setters ##
## No setters as users should not manipulate a 'nsmodelmoments' object ##
