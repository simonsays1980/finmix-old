#' @include datamoments.R
#' @include data.R
setClass("ddatamoments",
	representation(
		fact.moments = "matrix",
		over = "matrix",
		zeros = "matrix"
		),
	contains = c("datamoments"),
	prototype = list(
		fact.moments = matrix(), 
		over = matrix(), 
		zeros = matrix(), 
		data = data(), 
		mean = matrix(), 
		variance = matrix()),
	validity = function(object) {
		mom.moments <- object@fact.moments
		mom.data <- object@data
		mom.r <- getR(mom.data)
		if(mom.r != ncol(object@mean)) 
			return("Data dimension and dimension of the mean do not match.")
		if(mom.r != ncol(object@variance)) 
			return("Data dimension and dimension of the variance do not match.")
		if(mom.r != ncol(mom.moments) || nrow(mom.moments) != 4) 
			return("Data dimension and dimension of the L=4 factorial moments do not match.")
		if(mom.r != ncol(object@over)) 
			return("Data dimension and dimension of the overdispersion vector do not match.")
		if(mom.r != ncol(object@zeros)) 
			return("Data dimension and dimension of the zeros vector do not match.")
		if(!is.na(object@over) && any(object@over)) 
			return("Overdispersion is negative.")
		if(!is.na(object@zeros) && any(object@zeros < 0))
			return("Number of zeros is negative is negative.")
		if(!is.na(object@variance) && any(object@variance < 0))
			return("Variance is negative")
		## else: ok
		TRUE
	}
)

setMethod("initialize", "ddatamoments", function(.Object, ..., data) {
						.Object@data <- data
						dmoments(.Object) <- .Object@data
						callNextMethod(.Object, ...)
						return(.Object)
					}
)

"dmoments<-" <-  function(.Object, value) {
			              
				## Compute means ##
				## work only with data ordered by column ##
				if(!getByColumn(value)) {
					datam <- t(getY(value))
				}
				else {
					datam <- getY(value)
				}
				## means is a matrix 1 x r ##
				means <- matrix(0, ncol = getR(value), nrow = 1)
				for(i in 1:getR(value)) {
					means[i] <- mean(datam[,i])
				}
				.Object@mean <- means
				
				## Compute variances or variance-covariance matrix ##
				## variance is a r x r matrix ##
				.Object@variance <- matrix(diag(var(datam)), ncol = getR(value), nrow = 1)
	
				## Compute factorial moments ##
				## fact.moments is a L x r matrix (L = 4) ## 
				momentsm <- matrix(0, ncol = getR(value), nrow = 4)
				momentsm[1, ] <- means
				for(i in 1:getR(value)) {
					fact <- datam[,i] * max(datam[,i] - 1, 0)					
					momentsm[2, i] <- mean(fact)
					fact <- fact * max(datam[,i] - 2, 0)
					momentsm[3, i] <- mean(fact)
					fact <- fact * max(datam[,i] - 3, 0)
					momentsm[4, i] <- mean(fact)
				}
				.Object@fact.moments <- momentsm

				## Overdispersions and fractions of zeros ##
				## over and zeros are 1 x r matrices ## 
				zerosm <- matrix(0, ncol = getR(value), nrow = 1)
				.Object@over <- .Object@variance - as.matrix(means)
				for(i in 1:getR(value)) {
					zerosm[i] <- length(datam[datam[,i] == 0,i])/length(datam[,i])
				}
				.Object@zeros <- zerosm

				return(.Object)
}

setMethod("show", "ddatamoments", function(object) {
						.Object <- object
						dataname <- getName(getData(.Object))
						name <- ifelse(length(dataname) == 0, "", dataname)
						cat("Moments object of Data '", name, "'\n")
						colnames <- colnames(getY(getData(.Object)))
						if(!is.null(colnames)) {
							cat("			", colnames, "\n")
						}
						cat("	Mean			:[", getMean(.Object), "]\n")
						cat("	Variance		:[", getVariance(.Object), "]\n")
						cat("	Factorial Moments	:\n")
						hmom <- getFactorialMoments(.Object)
						for(i in 1:4) {
							cat("			 	[", hmom[i,], "]\n")
						}
						cat("	Overdispersion		:[", getOverDispersion(.Object), "]\n")
						cat("	Fractions of Zeros	:[", getZeros(.Object), "]\n")
						cat("	Data		:\n\n")
						cat("			 ", show(getData(.Object)), "\n")
					}
)

## Getters ##
setMethod("getMean", "ddatamoments", function(.Object) {
						return(.Object@mean)
					}
)
setMethod("getVariance", "ddatamoments", function(.Object) {
						return(.Object@variance)
					}
)
setMethod("getData", "ddatamoments", function(.Object) {
				return(.Object@data)
			}
)
setGeneric("getFactorialMoments", function(.Object) standardGeneric("getFactorialMoments"))
setMethod("getFactorialMoments", "ddatamoments", function(.Object) {
						return(.Object@fact.moments)
					}
)
setGeneric("getOverDispersion", function(.Object) standardGeneric("getOverDispersion"))
setMethod("getOverDispersion", "ddatamoments", function(.Object) {
						return(.Object@over)
					}
)
setGeneric("getZeros", function(.Object) standardGeneric("getZeros"))
setMethod("getZeros", "ddatamoments", function(.Object) {
						return(.Object@zeros)
					}
)
## Setters ##
setReplaceMethod("setMean", "ddatamoments", function(.Object, value) {
							.Object@mean <- value
							validObject(.Object)
							return(.Object)
						}
)

setReplaceMethod("setVariance", "ddatamoments", function(.Object, value) {
							.Object@variance <- value
							validObject(.Object)
							return(.Object)
						}
)
setReplaceMethod("setData", "ddatamoments", function(.Object, value) {
							.Object@data <- value
							moments(.Object) <- getData(.Object)
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setFactorialMoments<-", function(.Object, value) standardGeneric("setFactorialMoments<-"))			
setReplaceMethod("setFactorialMoments", "ddatamoments", function(.Object, value) {
							.Object@fact.moments <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setOverDispersion<-", function(.Object, value) standardGeneric("setOverDispersion<-"))			
setReplaceMethod("setOverDispersion", "ddatamoments", function(.Object, value) {
							.Object@over <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setZeros<-", function(.Object, value) standardGeneric("setZeros<-"))			
setReplaceMethod("setZeros", "ddatamoments", function(.Object, value) {
							.Object@zeros <- value
							validObject(.Object)
							return(.Object)
						}
)

