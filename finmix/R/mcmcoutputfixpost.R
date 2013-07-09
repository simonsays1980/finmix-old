setClass("mcmcoutputfixpost",
	representation(
		post 	= "list"),
	contains = c("mcmcoutputfix"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

setGeneric("getPost", function(object) standardGeneric("getPost"))
setMethod("getPost", "mcmcoutputfixpost", function(object) {
							return(object@post)	
						}
)


