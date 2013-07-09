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
	ST = "array",
	S = "array",
	NK = "array",
	post = "list",
	clust = "array",
	model = "model",
	prior = "prior"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

"mixturemcmc" <- function(data, model, prior, mcmc) {
	## check if all arguments are provided ##
	if(nargs() < 4) 
		return("[Error] All arguments must be provided.")

	## check data object ##
	validObject(data)
	has.data <- !all(is.na(data@y))
	has.S <- !all(is.na(data@S))
	K <- model@K
	if(has.data) {
		if(data@bycolumn) {
			datam <- data@y
			if(has.S && K == 1) {
				classm <- matrix(1, nrow = nrow(datam), ncol = 1)
			}
			else if(has.S && K > 1) {
				classm <- data@S
			}
		}
		else { ## data stored by row
			datam <- t(data@y)
			if(has.S && K == 1) {
				classm <- matrix(1, nrow = nrow(datam), ncol = 1)
			}
			else if(has.S && K > 1) {
				classm <- t(data@S)
			}
		}
		r <- ncol(datam)
		N <- nrow(datam)
	}
	else { ## data has no observations
		return("[Error] Observations in 'y' of 'data' object are obligatory for MCMC sampling.")
	}
	if(model@dist == "poisson") {
		has.exposures <- !all(is.na(data@exp))
		if(has.exposures) {
			if(data@bycolumn) {
				if(nrow(data@exp) != N && nrow(exposures) != 1) {
					return("[Error] number of exposures 'exp' in 'data' does not match number of observations in 'y'.")
				}
				else if(nrow(data@exp) == N) {
					exp <- data@exp
				}
				else { ## exp has dimension 1 x 1
					exp <-matrix(data@exp[1], nrow = N, ncol = 1)
				}
			}
			else { ## data stored by row
				if(ncol(data@exp) != N && ncol(data@exp) != 1) {
					return("[Error] number of exposures 'exp' in 'data' does not match number of observations in 'y'.")
				}
				else if(ncol(data@exp) == N) {
					exp <- t(data@exp)
				}
				else {
					exp <- matrix(data@exp[1], nrow = N, ncol = 1)
				}
			}
		}
		else { ## no exposures set default all to 1
			exp <- matrix(1, nrow = N, ncol = 1)
		}	
	}
	if(model@dist == "binomial") {
		## check for repetitions ##
		has.reps <- !all(is.na(data@T))
		if(has.reps) {
			if(bycolumn) {
				if(nrow(data@T) != N && nrow(data@T) != 1)  {
					return("[Error] number of repetitions 'T' in 'data' does not match number of observations in 'y'.")
				}
				else if(nrow(data@T) == N) {
					T <- data@T
				}
				else { ## dimension of T is 1 x 1
					T <- matrix(data@T[1], nrow = N, ncol = 1)
				}
			}
			else { ## data stored by row 
				if(ncol(data@T) != N && ncol(data@T) != 1) {
					return("[Error] number of repetitions 'T' in 'data' does not match number of observations in 'y'.")
				}
				else if(ncol(data@T) == N) {
					T <- t(data@T)
				}
				else { ## dimension of T is 1 x 1
					T <- matrix(data@T[1], nrow = N, ncol = 1)	
				}
			}
		}
		else { ## then check in model ##
			if(nrow(model@T) != N && nrow(model@T) != 1) {
				return("[Error] neither 'data' nor 'model' has correctly specified repetitions 'T'.")
			}	
			else if(nrow(model@T) == N) {
				T <- model@T
			}
			else { ## dimension of T is 1 x 1 
				T <- matrix(model@T[1], nrow = N, ncol = 1)
			}
		}
	}
	## check 'model' object ##
	validObject(model)
	K <- model@K
	norstud <- (model@dist == "normal" || model@dist == "normult" || model@dist == "student" || model@dist == "studmult")
	if(model@indicfix) {
		ranperm <- FALSE
	}
	if (mcmc@startpar && !model@indicfix && K > 1) { ## i.e. it should be started by sampling allocations
		if(length(model@par) == 0) 
			return("[Error] for starting with sampling allocations 'model' must provide starting parameters.")
		if(any(is.na(model@weight)) && model@indicmod == "multinomial") 
			return("[Error] for starting with sampling allocations 'model' must provide starting weights.")
	} 
	else { ## begin by sampling the parameters 
		if(!has.S && K > 1) 
			return("[Error] for starting with sampling parameters 'data' must provide starting allocations.")			
		if(norstud) {
			if(prior@type == "independent") { ## independent prior
				## later regression model ##

				## here only finite mixture ##
				if(length(model@par) == 0) 
					return("[Error] for an independent prior, starting values for the component means have to be provided.")
				has.mu <- "mu" %in% names(model@par)
				if(!has.mu) 
					return("[Error] for an independent prior, starting values for the component means have to be provided.")
			}
			norstudmult <- (model@dist == "normult" || model@dist == "studmult")
			has.logdet <- "logdet" %in% prior@par$sigma  
			if(norstudmult && !has.logdet) {
				has.C <- "C" %in% prior@par$sigma
				if(!has.C) { 
					return("[Error] for an independent prior entry 'C' in 'prior@par$sigma' has to be provided.")
				}
				else { ## if C is there check if array of dimension (r x r x K) 
					if(!is.array(C)) {
						return("[Error] 'C' in 'prior@par$sigma' must be an array of dimension (r x r x K).")
					}
					else {
						## check dimensions ##
						dims <- dim(prior@par$sigma$C)
						has.dim <- (dims[1] == r && dims[2] == r && dims[3] == K)
						if(!has.dim) 
							return("[Error] 'C' in 'prior@par$sigma' must be an array of dimension (r x r x K).")
						logdetC <- array(0, dim = c(r,r,K))
						for(k in 1:K) { ## TODO: check if Cs are given and check priordefine for thi
							logdetC[,,k] <- log(det(prior@par$sigma$C[,,k]))
						}
					}
				}
			} ## end norstudmult
		
		} ## end norstud 
			
		######################### MCMC SAMPLING #############################
		if(model@dist == "poisson") {
			M <- mcmc@M
			weights <- array(0, dim = c(M, K))
			pars <- list(lambda = array(0, dim = c(M, K)))
			if(prior@hier) {
				hypers <- list(b = array(0, dim = c(M, 1)))
			}
			else {
				hypers <- list(b = array(0, dim = c(1,1)))
			}
			log.mixlik <- array(0, dim = c(M, 1))
			log.mixprior <- array(0, dim = c(M, 1))
			if(!model@indicfix) {
				log.cdpost <- array(0, dim = c(M, 1))
				logs <- list(mixlik = log.mixlik, mixprior = log.mixprior, cdpost = log.cdpost)
			}
			else {
				logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
			}
			entropies <- array(0, dim = c(M, 1))
			if(!model@indicfix) {
				STm <- array(0, dim = c(M, 1))
				Sm <- array(integer(), dim = c(N, mcmc@storeS)) 
			} 
			else {
				STm <- array(0, dim = c(1, 1))
				Sm <- array(integer(), dim = c(1, 1))
			}
			NKs <- array(0, dim = c(M, K))
			if(mcmc@storepost) {
				postweights <- array(0, dim = c(M, K))
				posta <- array(0, dim = c(M, K))
				postb <- array(0, dim = c(M, K))
				postpars <- list(a = posta, b = postb)
				post <- list(par = postpars, weight = postweights)
			}
			else{
				postweights <- array(0, dim = c(1, K))
				posta <- array(0, dim = c(1, K))
				postb <- array(0, dim = c(1, K))
				postpars <- list(a = posta, b = postb)
				post <- list(par = postpars, weight = postweights)
			}
			if(!model@indicfix) {
				clustm <- array(0, dim = c(N, 1))
			}
			else {
				clustm <- array(0, dim = c(N, 1))
			}
			
			mcmcout <- new("mcmcoutput", name = data@name, M = mcmc@M, weight = weights, par = pars,
					ranperm = mcmc@ranperm, hyper = hypers, log = logs, entropy = entropies, 
					ST = STm, S = Sm, NK = NKs, post = post, model = model, prior = prior,
					clust = clustm)
			.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
			return(mcmcout)
		}	
	}
	
}
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
setGeneric("getST", function(.Object) standardGeneric("getST"))
setMethod("getST", "mcmcoutput", function(.Object) {
						return(.Object@ST)
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
setGeneric("getClust", function(.Object) standardGeneric("getClust"))
setMethod("getClust", "mcmcoutput", function(.Object) {
						return(.Object@clust)
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

