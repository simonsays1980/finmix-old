setClass("dataclass",
	representation(
	logpy = "matrix",
	prob = "matrix",
	mixlik = "numeric",
	entropy = "numeric",
	loglikcd = "matrix",
	postS = "numeric"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

"dataclass" <- function(data, model, simS = FALSE) {

		has.data <- !all(is.na(data@y))
		if(!has.data) {
			cat("[Error] dataclass needs observations 'y' in 'data' object to proceed.\n")
			return(FALSE)
		}
		if(model@indicfix) {
			cat("[Error] dataclass is not supposed to give back any values if allocations are given.\n")
			return(FALSE)
		}
		has.par <- !all(is.na(model@par))		
		if(!has.par) {
			cat("[Error] dataclass needs component parameters in 'model' object to proceed.\n")
			return(FALSE)
		}
		has.weight <- !all(is.na(model@weight)) 
		if (!has.weight) {
			cat("[Error] dataclass needs weight parameters in 'model' object to proceed.\n")
			return(FALSE)
		}
		## 'simS' for simulating S ##
		has.exposures <- !all(is.na(data@exp)) ## (exposures) can only be true for poisson
		has.T <- !all(is.na(data@T)) ## (repetitions) can only be true for binomial
		has.S <- !all(is.na(data@S))
		if(nargs() > 1) {
			if(data@bycolumn) {
				datam <- data@y
				N <- nrow(data@y)
				r <- data@r
				if(has.exposures) {
					exp <- data@exp
				}
				if(has.T) {
					Ti <- data@T
				}
				if(has.S) {
					classm <- data@S
				}
			}
			else {
				datam <- t(data@y)
				N <- nrow(data@y)
				if(has.exposures) {
					exp <- t(data@exp)
				}
				if(has.T) {
					Ti <- t(data@T)
				}
				if(has.S) {
					classm <- t(data@S)
				}
			}
		}	
		K <- model@K
		dist <- model@dist

		## Compute the likelihood l(y_i|theta_k) for all i and k ##		
		## lik.list is a 'list' object containing ##
                        ## 'lh'          exp(llh - maxl), an N x K  'matrix' 
                        ## 'maxl'        the maximum likelihood, an 1 x K 'vector' 
                        ## 'llh'         the likelihood, a N x K 'matrix' 
	
		## normal mixtures ##
		if(dist == "normal") {
			lik.list <- likelihood.normal(datam, model@par$mu, model@par$sigma)
		}

		## student-t mixtures ##
		else if(dist == "student") {
			lik.list <- likelihood.student(datam, model@par$mu, model@par$sigma, model@par$df)
		}
	
		## exponential mixtures ##
		else if(dist == "exponential") {
			lik.list <- likelihood.exponential(datam, model@par$lambda)
		}
	
		## poisson mixtures ##
		else if(dist == "poisson") {
			if(has.exposures) { ## should give a N x K 'matrix' object
				mu <- t(apply(exp, 1, "*", model@par$lambda))
			}	
			else { ## no exposures: 1 x K 'matrix' object
				mu <- model@par$lambda
			}
			lik.list <- likelihood.poisson(datam, mu)
		}
		
		## binomial mixtures ##
		else if(dist == "binomial") {
			lik.list <- likelihood.binomial(datam, Ti, model@par$p)					
		}

		## multivariate normal mixtures ##
		else if(dist == "normult") {
			has.sigmainv <- "sigmainv" %in% names(model@par)
			has.logdet <- "logdet" %in% names(model@par)
			if(has.sigmainv && has.logdet) {
				lik.list <- likelihood.normult(datam, model@par$mu, 
						model@par$sigmainv, model@par$logdet)
			}
			else {
				qinv <- array(0, dim = c(r,r,K))
				logdetq <- array(0, dim = c(1, K))
				for(k in 1:K) {
					qinv[,,k] <- solve(model@par$sigma[,,k])
					logdetq[k] <- log(det(qinv[,,k]))
				}
				lik.list <- likelihood.normult(datam, model@par$mu, 
						qinv, logdetq)
			}
		}
		
		## multivariate student-t mixtures ##
		else if(dist == "studmult") {
			has.sigmainv <- "sigmainv" %in% names(model@par)
			has.logdet <- "logdet" %in% names(model@par)
			if(has.sigmainv && has.logdet) {
				lik.list <- likelihood.studmult(datam, model@par$mu,
						model@par$sigmainv, model@par$logdet,
						model@par$df)
			}
			else {
				qinv <- array(0, dim = c(r,r,K))
				logdetq <- array(0, dim = c(1, K))
				for(k in 1:K) {
					qinv[,,k] <- solve(model@par$sigma[,,k])
					logdetq[k] <- log(det(qinv[,,k]))
				}
				lik.list <- likelihood.studmult(datam, model@par$mu,
						qinv, logdetq, model@par$df)
			}
		}

		## check attributes 'indicmod' and 'indicfix' in 'model' argument ##

		if(model@indicmod != "multinomial") {
			model@indicmod <- "multinomial"
		}	
		
		## Determine the posterior of the indicators S ##
		if(!model@indicfix) {
			if(model@indicmod == "multinomial") {
				
				## p is an N x K matrix ##
				p <- t(apply(lik.list$lh, 1, "*", model@weight))
				## sump is an N x 1 matrix ##
				sump <- apply(p, 1, sum)
				## lsump is an N x 1 matrix ##
				lsump <- log(sump) +  lik.list$maxl
				mixlik <- sum(lsump) ## numeric
				## p is the N x K probability classification matrix ##
				p = apply(p, 2, "/", sump)
				
				if(simS) {
					## simulate classifications from classification probability matrix ##
					if(K > 1) {
						rnd <- matrix(runif(nrow(datam)))
						S <- t(apply(p, 1, cumsum)) < matrix(rnd, nrow = nrow(rnd), ncol = K) 
						S <- matrix(apply(S, 1, sum)) + 1
						Sm <- matrix(S, nrow = nrow(S), ncol = K)
						Compm <- matrix(1:K, nrow = nrow(S), ncol = K, byrow = TRUE)
					        postS <- sum(log(apply((Sm == Compm) * p, 1, sum)))
					}
				}
	
				## compute complete data likelihood in case indicators were not simulated ##
				if(!simS) {
					if(has.S && K > 1) {
						loglikcd <- matrix(0, nrow = 1, ncol = K)
						for(k in 1:K) {
							loglikcd[k] <- sum(lik.list$llh[classm == k, k])
						}
					}
					else { ## no indicators given or no mixture ##
						loglikcd <- matrix(mixlik)
					}	
				}
			
			}
			## else: implemented later: Markov model for S ##
			
			## compute entropy ## 
			logp <- matrix(0, nrow = N, ncol = K)
			for(k in 1:K) {
				logp[p[,k] == 0,k] <- -99
				logp[p[,k] != 0,k] <- log(p[p[,k] != 0, k])
			}
			
			entropy <- -sum(logp * p)
			if(simS) {
				dataclass = new("dataclass", logpy = lik.list$llh, prob = p, mixlik = mixlik,
						entropy = entropy, loglikcd = matrix(), postS = postS)
				l <- list(dataclass = dataclass, S = as.integer(S))
				return(l)
			}
			dataclass <- new("dataclass", logpy = lik.list$llh, prob = p, mixlik = mixlik, 
			entropy = entropy, loglikcd = loglikcd, postS = numeric())

			return(dataclass)

		} 
		else { ## indicfix == TRUE
			if(has.S && K > 1) {
				loglikcd <- matrix(0, nrow = 1, ncol = K)
				for(k in 1:K) {
					loglikcd[k] <- sum(lik.list$llh[classm == k, k])
				}
			}
			else { ## no indicators given or no mixture ##
				loglikcd <- matrix(mixlik)
			}	

			dataclass = new("dataclass", logpy = lik.list$llh, prob = matrix(), mixlik = numeric(), 
						entropy = numeric(), loglikcd = loglikcd, postS = numeric())
			return(dataclass)
		}	
	
}

setMethod("show", "dataclass", function(object) {
					has.loglikcd <- !all(is.na(object@loglikcd))
					has.postS <- !all(is.na(object@postS))
					cat("Dataclass object\n")
					cat("	logpy		:", paste(dim(object@logpy), collapse = "x"), "\n")
					cat("	prob		:", paste(dim(object@prob), collapse = "x"), "\n")
					cat("	mixlik		:", object@mixlik, "\n")
					cat("	entropy		:", object@entropy, "\n")
					if(has.loglikcd) {
						cat("	loglikcd	:", paste(dim(object@loglikcd), collapse = "x"), "\n")	
					}
					if(has.postS) {
						cat("	postS		:", object@postS, "\n")
					}
				}
)
## Getters ##
setGeneric("getLogPY", function(.Object) standardGeneric("getLogPY"))
setMethod("getLogPY", "dataclass", function(.Object) {
						return(.Object@logpy)
					}
)
setGeneric("getProb", function(.Object) standardGeneric("getProb"))
setMethod("getProb", "dataclass", function(.Object) {
						return(.Object@prob)
					}
)
setGeneric("getMixlik", function(.Object) standardGeneric("getMixlik"))
setMethod("getMixlik", "dataclass", function(.Object) {
						return(.Object@mixlik)
					}
)
setGeneric("getEntropy", function(.Object) standardGeneric("getEntropy"))
setMethod("getEntropy", "dataclass", function(.Object) {
						return(.Object@entropy)
					}
)
setGeneric("getLoglikcd", function(.Object) standardGeneric("getLoglikcd"))
setMethod("getLoglikcd", "dataclass", function(.Object) {
						return(.Object@loglikcd)
					}
)
setGeneric("getPostS", function(.Object) standardGeneric("getPostS"))
setMethod("getPostS", "dataclass", function(.Object) {
						return(.Object@postS)
					}
)

## Setters ##
setGeneric("setLogPY<-", function(.Object, value) standardGeneric("setLogPY<-"))
setReplaceMethod("setLogPY", "dataclass", function(.Object, value) {
						.Object@logpy <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setProb<-", function(.Object, value) standardGeneric("setProb<-"))
setReplaceMethod("setProb", "dataclass", function(.Object, value) {
							.Object@prob <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setMixlik<-", function(.Object, value) standardGeneric("setMixlik<-"))
setReplaceMethod("setMixlik", "dataclass", function(.Object, value) {
							.Object@mixlik <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setEntropy<-", function(.Object, value) standardGeneric("setEntropy<-"))
setReplaceMethod("setEntropy", "dataclass", function(.Object, value) {
							.Object@entropy <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setLoglikcd<-", function(.Object, value) standardGeneric("setLoglikcd<-"))
setReplaceMethod("setLoglikcd", "dataclass", function(.Object, value) {
							.Object@loglikcd <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setPostS<-", function(.Object, value) standardGeneric("setPostS<-"))
setReplaceMethod("setPostS", "dataclass", function(.Object) {
							.Object@postS <- value
							validObject(.Object)
							return(.Object)
						}
)

