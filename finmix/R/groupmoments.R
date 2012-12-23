#' @include data.R
setClass("groupmoments", 
	representation(
		NK = "matrix",
		group.mean = "matrix",
		g.data = "data"), 
	prototype = list(
		NK = matrix(),
		group.mean = matrix(),
		g.data = data()
		),
	validity = function(object) {
				mom.nk <- object@NK
				mom.group.means <- object@group.mean
				mom.K <- length(levels(factor(getY(object@g.data))))
				if(mom.K != ncol(mom.nk)) 
					return("Data dimension and dimension of group sizes do not match.")
				if(mom.K != ncol(mom.group.means))
					return("Data dimension and dimension of the group means fo not match.")
				if(any(mom.nk < 0))
					return("Group size is negative.") 
				## else: ok
				TRUE
			}
)

setMethod("initialize", "groupmoments", function(.Object, ..., data) {
							.Object@g.data <- data
							gmoments(.Object) <- .Object@g.data
							callNextMethod(.Object, ...)
							return(.Object)
						}
)

"gmoments<-" <- function(.Object, value) {
			## Compute group sizes ##
			## work only with  ordered by column ##
			if(!getByColumn(value)) {
				datam <- t(getY(value))
				classm <- t(getS(value))
			}
			else {
				datam <- getY(value)
				classm <- getS(value)
			}
			
			## Calculate group sizes and group means ##
			## 'NK' is a 1 x K matrix ##
			## 'group.mean' is a 1 x K matrix for r == 1 ## 
			## 'group.mean' is a r x K matrix for r > 1 ##
			level.set <- as.numeric(levels(factor(classm)))
			if(getR(value) == 1) {
				nkm <- matrix(0, ncol = length(level.set), nrow = 1)
				group.means <- matrix(0, ncol = length(level.set), nrow = 1)
				for(i in 1:length(level.set)) {
					nkm[i] <- sum(match(classm, level.set[i]), na.rm = TRUE)
					group.means[i] <- mean(datam[which(classm == level.set[i]), drop = FALSE], na.rm = TRUE)
				}
			}
			else {
				nkm <- matrix(0, ncol = length(level.set), nrow = 1)
				group.means <- matrix(0, ncol = length(level.set), nrow = getR(value))
				for(i in 1:length(level.set)) {
					nkm[i] <- sum(match(classm, level.set[i]), na.rm = TRUE)
					group.means[,i] <- apply(datam[which(classm == level.set[i]), , drop = FALSE], 2, mean,na.rm = TRUE)
				}
			}
			.Object@NK <- as.matrix(nkm)
			.Object@group.mean <- as.matrix(group.means)

			
			return(.Object)
}

## R usual 'show' function ##
setMethod("show", "groupmoments", function(object) {
						data.name <- getName(getG.Data(object))
						oname <- ifelse(length(data.name) == 0, "", data.name)
						cat("Groupmoments object '", oname, "'\n")
						cat("	NK		: [", format(getNK(object), trim = TRUE), "]\n")
						cat("	group.mean	: [", format(getGroupMean(object)[1,], trim = TRUE), "]\n")
						r <- nrow(getGroupMean(object))
						if(r > 1) {
							for(i in 2:r) 
								cat("		  [", format(getGroupMean(object)[i,], trim = TRUE), "]\n")
						}
						
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
setGeneric("getG.Data", function(.Object) standardGeneric("getG.Data"))
setMethod("getG.Data", "groupmoments", function(.Object) {
						return(.Object@g.data)		
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
setGeneric("setG.Data<-", function(.Object, value) standardGeneric("setG.Data<-"))
setReplaceMethod("setG.Data", "groupmoments", function(.Object, value) {
							.Object@g.data <- value
							validObject(.Object)
							return(.Object)
						}
)
