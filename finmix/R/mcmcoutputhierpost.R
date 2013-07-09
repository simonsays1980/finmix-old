setClass("mcmcoutputhierpost", 
	contains = c("mcmcoutputhier", "mcmcoutputpost"),
	validity = function(object) {
			## else: OK
			TRUE
	}
)

## Set 'mcmcoutput' to the virtual class inheriting 	##
## to each other 'mcmcoutput' class. 			##
## This is done to simplify dispatching methods.	##
setClassUnion("mcmcoutput", 
	c(
		"mcmcoutputfix",
		"mcmcoutputfixhier",
		"mcmcoutputfixpost",
		"mcmcoutputfixhierpost",
		"mcmcoutputbase",
		"mcmcoutputhier",
		"mcmcoutputpost",
		"mcmcoutputhierpost")
)
