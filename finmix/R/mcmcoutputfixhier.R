setClass("mcmcoutputfixhier",
	representation(
		hyper = "list"),
	contains = c("mcmcoutputfix"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)


setGeneric("getHyper", function(object) standardGeneric("getHyper"))
setMethod("getHyper", "mcmcoutputfixhier", function(object) {
							return(object@hyper)
						}
)
	
## No setters for this object as it is not intended 	##
## that users manipulate this object 		    	##
