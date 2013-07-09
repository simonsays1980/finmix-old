setClass("mcmcoutputhier",
	representation(
		hyper = "list"),
	contains = c("mcmcoutputbase"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

## Generic method already set in class 'mcmcoutputfixhier' ##
setMethod("getHyper", "mcmcoutputhier", function(object) {
					return(object@hyper)
				}
)
