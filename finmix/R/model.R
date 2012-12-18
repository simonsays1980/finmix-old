setClass("model",
	representation(
		K="numeric",
		dist="character"),
	validity = function(object) {
				if( !is.numeric(object@K) ) 
					return("Number of mixture components 'K' has to be of class numeric!")
				## check if the distribution is part of a set of implemented distribution functions ##
				## else: OK ##
				TRUE
			},
	prototype = list(K = 2, dist = "Normal")

)

## Constructor for class 'model' ##
"model" <- function(K=2, dist="Normal") {
		model <- new("model", K = K, dist = dist)
		return(model)
}

  
