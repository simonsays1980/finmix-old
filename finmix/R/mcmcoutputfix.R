setClass("mcmcoutputfix",
	representation(
		M 	= "integer",
		ranperm = "logical",
		par 	= "list",
		log	= "list",
		model 	= "model",
		prior	= "prior"),
	validity = function(object) {
			##else: OK
			TRUE
	}
)

setGeneric("getM", function(object) standardGeneric("getM")) 
setMethod("getM", "mcmcoutputfix", function(object) {
						return(object@M)
					}
)
setGeneric("getRanPerm", function(object) standardGeneric("getRanPerm"))
setMethod("getRanPerm", "mcmcoutputfix", function(object) {
							return(object@ranperm)
						}
)
## Generic set in model.R ##
setMethod("getPar", "mcmcoutputfix", function(.Object) {
						return(.Object@par)
					}
)
setGeneric("getLog", function(object) standardGeneric("getLog"))
setMethod("getLog", "mcmcoutputfix", function(object) {
						return(object@log)
					}
)

setGeneric("getModel", function(object) standardGeneric("getModel"))
setMethod("getModel", "mcmcoutputfix", function(object) {
						return(object@model)
					}
)
setGeneric("getPrior", function(object) standardGeneric("getPrior"))
setMethod("getPrior", "mcmcoutputfix", function(object) {
						return(object@prior)
					}
)

## no setters: users should not get access to manipulate ##
## data from MCMC output				 ##
