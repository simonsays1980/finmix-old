setClass("modelmoments", 
	representation(
	mean = "matrix",
	var = "matrix"
	),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

"modelmoments" <- function(model, J. = 4) {
		dist <- model@dist
		if(dist == "normal" || dist == "normult" ) {
			
			nsmodelmoments <- nsmodelmoments(model, J.)
			return(nsmodelmoments)
		}
		else if(dist == "student" || dist == "studmult") {
			return("[Warning] Function 'modelmoments' not implemented for Student mixtures.")
		}
		else if(dist == "poisson") {
			
			dmodelmoments <- dmodelmoments(model, J.)
			return(dmodelmoments)
		}
		else if(dist == "binomial") {
			return("[Warning] Function 'modelmoments' not implemented for binomial mixtures.")
		}
		else if(dist == "exponential") {
			return("[Warning] Function 'modelmoments' not implemented for exponential mixtures")
		}
		
}

## Getters ##
setGeneric("getMean", function(.Object) standardGeneric("getMean"))
setMethod("getMean", "modelmoments", function(.Object) {
						return(.Object@mean)
					}		
)
setGeneric("getVar", function(.Object) standardGeneric("getVar"))
setMethod("getVar", "modelmoments", function(.Object) {
						return(.Object@var)
					}
)

## Setters are not provided as users should not manipulate a 'modelmoments' object ##
