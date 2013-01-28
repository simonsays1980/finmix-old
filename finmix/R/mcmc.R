setClass("mcmc", 
	representation(
	burnin = "integer",
	M = "integer",
	startpar = "logical",
	storeS = "integer",
	storepost = "logical",
	ranperm = "logical"),
	validity = function(object) {
		if(object@burnin < 0)
			return("[Error] Number of Burn-In draws 'burnin' is negative.")
		if(object@M < 0)
			return("[Error] Number of draws 'M' is negative.")
		if(object@M == 0)
			return("[Error] Number of draws 'M' is zero.")
		if(object@storeS < 0)
			return("[Error] Number of draws of S to store 'storeS' is negative.")
		## else: OK
		TRUE
	}
)	
"mcmc" <- function(burnin. = 0, M. = 5000, startpar. = FALSE, storeS. = 1000, storepost. = TRUE,
			ranperm. = TRUE) {
		
		mcmc <- new("mcmc", burnin = burnin., M = M., startpar = startpar., storeS = storeS.,
				storepost = storepost., ranperm = ranperm.)

		return(mcmc)
}

## Getters ##
setGeneric("getBurnin", function(.Object) standardGeneric("getBurnin"))
setMethod("getBurnin", "mcmc", function(.Object) {
					return(.Object@burnin)
				}
)
setGeneric("getM", function(.Object) standardGeneric("getM"))
setMethod("getM", "mcmc", function(.Object) {
					return(.Object@M)
				}
)
setGeneric("getStartpar", function(.Object) standardGeneric("getStartpar"))
setMethod("getStartpar", "mcmc", function(.Object) {
					return(.Object@startpar)
				}
)
setGeneric("getStoreS", function(.Object) standardGeneric("getStoreS"))
setMethod("getStoreS", "mcmc", function(.Object) {
					return(.Object@storeS)
				}
)
setGeneric("getStorepost", function(.Object) standardGeneric("getStorepost"))
setMethod("getStorepost", "mcmc", function(.Object) {
					return(.Object@storepost)
				}
)
setGeneric("getRanperm", function(.Object) standardGeneric("getRanperm"))
setMethod("getRanperm", "mcmc", function(.Object) {
					return(.Object@ranperm)
				}
)

## Setters ##
setGeneric("setBurnin<-", function(.Object, value) standardGeneric("setBurnin<-"))
setReplaceMethod("setBurnin", "mcmc", function(.Object, value) {
						.Object@burnin <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setM<-", function(.Object, value) standardGeneric("setM<-"))
setReplaceMethod("setM", "mcmc", function(.Object, value) {
						.Object@M <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setStartpar<-", function(.Object, value) standardGeneric("setStartpar<-"))
setReplaceMethod("setStartpar", "mcmc", function(.Object, value) {
						.Object@startpar <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setStoreS<-", function(.Object, value) standardGeneric("setStoreS<-"))
setReplaceMethod("setStoreS", "mcmc", function(.Object, value) {
						.Object@storeS <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setStorepost<-", function(.Object, value) standardGeneric("setStorepost<-"))
setReplaceMethod("setStorepost", "mcmc", function(.Object, value) {
						.Object@storepost <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setRanperm<-", function(.Object, value) standardGeneric("setRanperm<-"))
setReplaceMethod("setRanperm", "mcmc", function(.Object, value) {
						.Object@ranperm <- value
						validObject(.Object)
						return(.Object)
					}
)

