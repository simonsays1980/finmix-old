setClass("mcmcoutputbase",
	representation(
		weight 	= "array",
		entropy	= "array",
		ST 	= "array",
		S 	= "array",
		NK 	= "array",
		clust 	= "array"),
	contains = c("mcmcoutputfix"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

setGeneric("getWeight", function(object) standardGeneric("getWeight"))
setMethod("getWeight", "mcmcoutputbase", function(object) {
							return(object@weight)
						}
)
setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))
setMethod("getEntropy", "mcmcoutputbase", function(object) {
							return(object@entropy)	
						}
)
setGeneric("getST", function(object) standardGeneric("getST"))
setMethod("getST", "mcmcoutputbase", function(object) {
							return(object@ST)	
						}
)
setGeneric("getS", function(object) standardGeneric("getS"))
setMethod("getS", "mcmcoutputbase", function(object) {
							return(object@S)	
						}
)
setGeneric("getNK", function(object) standardGeneric("getNK"))
setMethod("getNK", "mcmcoutputbase", function(object) {
							return(object@NK)	
						}
)
setGeneric("getClust", function(object) standardGeneric("getClust"))
setMethod("getClust", "mcmcoutputbase", function(object) {
							return(object@clust)	
						}
)


