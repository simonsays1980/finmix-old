setClass("mcmcoutputpost",
	representation(
		post = "list"),
	contains = c("mcmcoutputbase"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

setMethod("getPost", "mcmcoutputpost", function(object) {
						return(object@post)
					}
)
