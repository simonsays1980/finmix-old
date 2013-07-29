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

setClass("groupmoments", 
	representation(
		NK = "matrix",
		group.mean = "matrix",
		data = "data"), 
	validity = function(object) {
				mom.nk <- object@NK
				mom.group.means <- object@group.mean
				mom.K <- length(levels(factor(getS(object@data))))
				if(mom.K != ncol(mom.nk)) 
					return("[Error] Data dimension and dimension of group sizes do not match.")
				if(mom.K != ncol(mom.group.means))
					return("[Error] Data dimension and dimension of the group means fo not match.")
				if(any(mom.nk < 0))
					return("[Error] Group size is negative.") 
				## else: ok
				TRUE
			}
)

"groupmoments" <- function(data) {
			## Compute group sizes ##
			## work only with  ordered by column ##
			if(!data@bycolumn) {
				datam <- t(data@y)
				classm <- t(data@S)
			}
			else {
				datam <- data@y
				classm <- data@S
			}
			
			## Calculate group sizes and group means ##
			## 'NK' is a 1 x K matrix ##
			## 'group.mean' is a 1 x K matrix for r == 1 ## 
			## 'group.mean' is a r x K matrix for r > 1 ##
			level.set <- as.numeric(levels(factor(classm)))
			if(data@r == 1) {
				nkm <- matrix(0, ncol = length(level.set), nrow = 1)
				group.means <- matrix(0, ncol = length(level.set), nrow = 1)
				for(i in 1:length(level.set)) {
					nkm[i] <- sum(match(classm, level.set[i]), na.rm = TRUE)
					group.means[i] <- mean(datam[which(classm == level.set[i]), drop = FALSE], na.rm = TRUE)
				}
			}
			else {
				nkm <- matrix(0, ncol = length(level.set), nrow = 1)
				group.means <- matrix(0, ncol = length(level.set), nrow = data@r)
				for(i in 1:length(level.set)) {
					nkm[i] <- sum(match(classm, level.set[i]), na.rm = TRUE)
					group.means[,i] <- apply(datam[which(classm == level.set[i]), , drop = FALSE], 2, mean,na.rm = TRUE)
				}
			}
		
			groupmoments <- new("groupmoments", NK = as.matrix(nkm), group.mean = as.matrix(group.means), 
						data = data)
			
			return(groupmoments)
}

## R usual 'show' function ##
setMethod("show", "groupmoments", function(object) {
						data.name <- getName(object@data)
						oname <- ifelse(length(data.name) == 0, "", data.name)
						cat("Groupmoments object '", oname, "'\n")
						cat("	Type		:", class(object), "\n")
						cat("	NK		:", paste(dim(object@NK), collapse="x"), "\n")
						cat("	Group Mean	:", paste(dim(object@group.mean), collapse="x"), "\n")
					##	cat("	NK		: [", format(object@NK, trim = TRUE), "]\n")
					##	cat("	group.mean	: [", format(object@group.mean[1,], trim = TRUE), "]\n")
					##	r <- nrow(object@group.mean)
					##	if(r > 1) {
					##		for(i in 2:r) 
					##			cat("		  	  [", format(object@group.mean[i,], trim = TRUE), "]\n")
					##	}
						
					}
)
 
## R usual Getters ##
setGeneric("getNK", function(.Object) standardGeneric("getNK"))
setMethod("getNK", "groupmoments", function(.Object) {
						return(.Object@NK)					
					}
)
setGeneric("getGroupMean", function(.Object) standardGeneric("getGroupMean"))
setMethod("getGroupMean", "groupmoments", function(.Object) {
						return(.Object@group.mean)	
					}
)
setGeneric("getData", function(.Object) standardGeneric("getData"))
setMethod("getData", "groupmoments", function(.Object) {
						return(.Object@data)		
					}
)
## R usual Setters ##
setGeneric("setNK<-", function(.Object, value) standardGeneric("setNK<-"))
setReplaceMethod("setNK", "groupmoments", function(.Object, value) {
							.Object@NK <- value
							return(.Object)							
						}
)
setGeneric("setGroupMean<-", function(.Object, value) standardGeneric("setGroupMean<-"))
setReplaceMethod("setGroupMean", "groupmoments", function(.Object, value) {
						.Object@group.mean <- value
						return(.Object)
					}
)
setGeneric("setData<-", function(.Object, value) standardGeneric("setData<-"))
setReplaceMethod("setData", "groupmoments", function(.Object, value) {
							.Object@data <- value
							validObject(.Object)
							return(.Object)
						}
)
