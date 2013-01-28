setClass("dmodelmoments", 
	representation(
	over = "numeric",
	factorial = "matrix",
	zero = "numeric"),
	contains = c("modelmoments"),	
	validity = function(object) {
		## else: OK
		TRUE
	}
)

"dmodelmoments" <- function(model, J.) {
		
			K <- model@K
			
			meanm <- matrix(sum(model@par$lambda * model@weight))
			varm <- matrix(sum(model@weight * model@par$lambda * (model@par$lambda + 1)) - meanm^2)
			if(K > 1) {
				over <- varm[1] - meanm[1]
			} 
			else {
				over <- 0
			}
			factorial <- matrix(0, ncol = J., nrow = 1)
				
			for(i in 1:J.) {
				factorial[i] <- sum(model@weight * model@par$lambda^i)
			}
			zero <- sum(model@weight * exp(-model@par$lambda))
			
			dmodemoments <- new("dmodelmoments", mean = meanm, var = varm, over = over, 
					factorial = factorial, zero = zero)
}

setMethod("show", "dmodelmoments", function(object) {
						cat("Modelmoments object\n")
						cat("	Type		:", class(object), "\n")
						cat("	Mean		:", object@mean, "\n")
						cat("	Var		:", object@var, "\n")
						cat("	Over		:", object@over, "\n")
						cat("	Factorial	:", paste(dim(object@factorial), collapse = "x"), "\n")
						cat("	Zero		:", object@zero, "\n")
					}
)
## Getters ##
## Generic set in 'modelmoments' class ##
setMethod("getMean", "dmodelmoments", function(.Object) {
						return(.Object@mean)
					}
)
## Generic set in 'modelmoments' class ##
setMethod("getVar", "dmodelmoments", function(.Object) {
						return(.Object@var)
					}
)
setGeneric("getOver", function(.Object) standardGeneric("getOver"))
setMethod("getOver", "dmodelmoments", function(.Object) {
						return(.Object@over)
					}
)
setGeneric("getFactorial", function(.Object) standardGeneric("getFactorial"))
setMethod("getFactorial", "dmodelmoments", function(.Object) {
							return(.Object@factorial)
						}
)
setGeneric("getZero", function(.Object) standardGeneric("getZero"))
setMethod("getZero", "dmodelmoments", function(.Object) {
						return(.Object@zero)
					}
)

## Setters ##
## No setters as users should not manipulate a 'dmodelmoments' object ##
