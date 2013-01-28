setClass("mcmcoutput",
	representation(
	name = "character",
	M = "integer",
	weight = "array",
	par = "list",
	ranperm = "logical",
	hyper = "list",
	log = "list",
	entropy = "array",
	S = "array",
	NK = "array",
	post = "list",
	model = "model",
	prior = "prior"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

## Getters ##
## Generic set in 'data' class ##
setMethod("getName", "mcmcoutput", function(.Object) {
						return(.Object@name)
					}
)
## Generic set in 'mcmc' class ##
setMethod("getM", "mcmcoutput", function(.Object) {
					return(.Object@M)
				}
)
## Generic set in 'model' class ##
setMethod("getWeight", "mcmcoutput", function(.Object) {
						return(.Object@weight)
					}
)
## Generic set in 'model' class ##
setMethod("getPar", "mcmcoutput", function(.Object) {
						return(.Object@par)
					}
)
## Generic set in 'mcmc' class ##
setMethod("getRanperm", "mcmcoutput", function(.Object) {
						return(.Object@ranperm)
					}
)
setGeneric("getHyper", function(.Object) standardGeneric("getHyper"))
setMethod("getHyper", "mcmcoutput", function(.Object) {
						return(.Object@hyper)
					}
)
setGeneric("getLog", function(.Object) standardGeneric("getLog"))
setMethod("getLog", "mcmcoutput", function(.Object) {
						return(.Object@log)
					}
)
setGeneric("getEntropy", function(.Object) standardGeneric("getEntropy"))
setMethod("getEntropy", "mcmcoutput", function(.Object) {
						return(.Object@entropy)
					}
)
## Generic set in 'data' class ##
setMethod("getS", "mcmcoutput", function(.Object) {
					return(.Object@S)
				}
)
## Generic set in 'groupmoments' class ##
setMethod("getNK", "mcmcoutput", function(.Object) {
						return(.Object@NK)
					}
)
setGeneric("getPost", function(.Object) standardGeneric("getPost"))
setMethod("getPost", "mcmcoutput", function(.Object) {
						return(.Object@post)
					}
)
setGeneric("getModel", function(.Object) standardGeneric("getModel"))
setMethod("getModel", "mcmcoutput", function(.Object) {
						return(.Object@model)		
					}
)
setGeneric("getPrior", function(.Object) standardGeneric("getPrior"))
setMethod("getPrior", "mcmcoutput", function(.Object) {
						return(.Object@prior)
					}
)

## Setters ##
## Only 'name' attribute is accessible by user as 'mcmcoutput' is an output object with results ##
## Generic already set in 'data' class ## 
setReplaceMethod("setName", "mcmcoutput", function(.Object, value) {
							.Object@name <- value
							validObject(.Object)
							return(.Object)
						}
)

