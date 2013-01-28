setClass("ddatamoments",
	representation(
		factorial = "matrix",
		over = "matrix",
		zero = "matrix"
		),
	contains = c("datamoments"),
	validity = function(object) {
		mom.moments <- object@factorial
		mom.data <- object@data
		mom.r <- getR(mom.data)
		if(mom.r != ncol(object@mean)) 
			return("[Error] Data dimension and dimension of the mean do not match.")
		if(mom.r != ncol(object@variance)) 
			return("[Error] Data dimension and dimension of the variance do not match.")
		if(mom.r != ncol(mom.moments) || nrow(mom.moments) != 4) 
			return("[Error] Data dimension and dimension of the L=4 factorial moments do not match.")
		if(mom.r != ncol(object@over)) 
			return("[Error] Data dimension and dimension of the overdispersion vector do not match.")
		if(mom.r != ncol(object@zero)) 
			return("[Error] Data dimension and dimension of the zeros vector do not match.")
		if(!is.na(object@over) && any(object@over)) 
			return("[Error] Overdispersion is negative.")
		if(!is.na(object@zero) && any(object@zero < 0))
			return("[Error] Number of zeros is negative is negative.")
		if(!is.na(object@variance) && any(object@variance < 0))
			return("[Error] Variance is negative")
		## else: ok
		TRUE
	}
)

setMethod("initialize", "ddatamoments", function(.Object, ..., data) {
						.Object@data <- data
						dmoments(.Object) <- .Object@data
						has.S <- !all(is.na(data@S))
						if(has.S) {
							.Object@sdatamoments <- sdatamoments(data)
						}
						callNextMethod(.Object, ...)
						return(.Object)
					}
)

"dmoments<-" <-  function(.Object, value) {
			              
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
				.Object@var <- matrix(diag(var(datam)), ncol = value@r, nrow = 1)
				if(has.colnames) {
					colnames(.Object@var) <- colnames(datam)
					rownames(.Object@var) <- colnames(datam)
				}

				## Compute factorial moments ##
				## fact.moments is a L x r matrix (L = 4) ## 
				momentsm <- matrix(0, nrow = value@r, ncol = 4)
				momentsm[1, ] <- means
				for(i in 1:value@r) {
					fact <- datam[,i] * max(datam[,i] - 1, 0)					
					momentsm[i, 2] <- mean(fact)
					fact <- fact * max(datam[,i] - 2, 0)
					momentsm[i, 3] <- mean(fact)
					fact <- fact * max(datam[,i] - 3, 0)
					momentsm[i, 4] <- mean(fact)
				}
				if(has.colnames) {
					rownames(momentsm) <- colnames(datam)
				}
				.Object@factorial <- momentsm

				## Overdispersions and fractions of zeros ##
				## over and zeros are r x 1 matrices ## 
				zerosm <- matrix(0, nrow = value@r, ncol = 1)
				.Object@over <- t(diag(.Object@var) - as.matrix(means))
				for(i in 1:value@r) {
					zerosm[i] <- length(datam[datam[,i] == 0,i])/length(datam[,i])
				}
				.Object@zero <- zerosm
				if(has.colnames) {
					rownames(.Object@over) <- colnames(datam)
					rownames(.Object@zero) <- colnames(datam)
				}

				return(.Object)
}

setMethod("show", "ddatamoments", function(object) {
						dataname <- getName(object@data)
						name <- ifelse(length(dataname) == 0, "", dataname)
						cat("Datamoments object '", name, "'\n")
						cat("	Type		:", class(object), "\n")
						cat("	Mean		:", paste(dim(object@mean), collapse = "x") , "\n")
						cat("	Var		:", paste(dim(object@var), collapse = "x"), "\n")
						cat("	Factorial	:", paste(dim(object@factorial), collapse = "x"), "\n")
						cat("	Over		:", paste(dim(object@over), collapse = "x"), "\n")
						cat("	Zero		:", paste(dim(object@zero), collapse = "x"), "\n")
						has.S <- !all(is.na(getS(getData(object))))
						if(has.S) {
							cat("	SDataMoments	:", class(object@sdatamoments), "\n")
						}
					}
)

## Getters ##
setMethod("getMean", "ddatamoments", function(.Object) {
						return(.Object@mean)
					}
)
setMethod("getVar", "ddatamoments", function(.Object) {
						return(.Object@var)
					}
)
setMethod("getData", "ddatamoments", function(.Object) {
				return(.Object@data)
			}
)
setMethod("getSDataMoments", "ddatamoments", function(.Object) {
					return(.Object@sdatamoments)
				}
)
setGeneric("getFactorial", function(.Object) standardGeneric("getFactorial"))
setMethod("getFactorial", "ddatamoments", function(.Object) {
						return(.Object@factorial)
					}
)
setGeneric("getOver", function(.Object) standardGeneric("getOver"))
setMethod("getOver", "ddatamoments", function(.Object) {
						return(.Object@over)
					}
)
setGeneric("getZero", function(.Object) standardGeneric("getZero"))
setMethod("getZero", "ddatamoments", function(.Object) {
						return(.Object@zero)
					}
)

## Setters ##
## No setters as users should not manipulate a 'ddatamoments' object ## 
