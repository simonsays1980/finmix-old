## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

"mixturemcmc" <- function(data, model, prior, mcmc) {
	## check if all arguments are provided ##
	if(nargs() < 4) 
		stop("All arguments must be provided.")

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
		stop("Observations in 'y' of 'data' object are obligatory for MCMC sampling.")
	}
	if(model@dist == "poisson" || model@dist == "cond.poisson") {
		has.exposures <- !all(is.na(data@exp))
		if(has.exposures) {
			if(data@bycolumn) {
				if(nrow(data@exp) != N && nrow(data@exp) != 1) {
					stop("Number of exposures 'exp' in 'data' does not match number of observations in 'y'.")
				}
				else if(nrow(data@exp) == N) {
					data@exp <- data@exp
				}
				else { ## exp has dimension 1 x 1
					data@exp <- matrix(data@exp[1, 1], nrow = N, ncol = 1)
				}
			}
			else { ## data stored by row
				if(ncol(data@exp) != N && ncol(data@exp) != 1) {
					stop("Number of exposures 'exp' in 'data' does not match number of observations in 'y'.")
				}
				else if(ncol(data@exp) == N) {
					data@exp <- t(data@exp)
				}
				else {
					data@exp <- matrix(data@exp[1, 1], nrow = N, ncol = 1)
				}
			}
		}
		else { ## no exposures set default all to 1
			data@exp <- matrix(1, nrow = N, ncol = 1)
		}	
	}
	if(model@dist == "binomial") {
		## check for repetitions ##
		has.reps <- !all(is.na(data@T))
		if(has.reps) {
			if(bycolumn) {
				if(nrow(data@T) != N && nrow(data@T) != 1)  {
					stop("Number of repetitions 'T' in 'data' does not match number of observations in 'y'.")
				}
				else if(nrow(data@T) == N) {
					T <- data@T
				}
				else { ## dimension of T is 1 x 1
					T <- matrix(data@T[1, 1], nrow = N, ncol = 1)
				}
			}
			else { ## data stored by row 
				if(ncol(data@T) != N && ncol(data@T) != 1) {
					stop("Number of repetitions 'T' in 'data' does not match number of observations in 'y'.")
				}
				else if(ncol(data@T) == N) {
					T <- t(data@T)
				}
				else { ## dimension of T is 1 x 1
					T <- matrix(data@T[1, 1], nrow = N, ncol = 1)	
				}
			}
		}
		else { ## then check in model ##
			if(nrow(model@T) != N && nrow(model@T) != 1) {
				stop("Neither 'data' nor 'model' has correctly specified repetitions 'T'.")
			}	
			else if(nrow(model@T) == N) {
				T <- model@T
			}
			else { ## dimension of T is 1 x 1 
				T <- matrix(model@T[1, 1], nrow = N, ncol = 1)
			}
		}
	}
	## check 'model' object ##
	validObject(model)
	K <- model@K
	norstud <- (model@dist == "normal" || model@dist == "normult" || model@dist == "student" || model@dist == "studmult")
	if (model@indicfix) {
		ranperm <- FALSE
	}
	if (mcmc@startpar && !model@indicfix && K > 1) { ## i.e. it should be started by sampling allocations
		if(length(model@par) == 0) 
			stop("For starting with sampling allocations 'model' must provide starting parameters.")
		if(any(is.na(model@weight)) && model@indicmod == "multinomial") 
			stop("For starting with sampling allocations 'model' must provide starting weights.")
		dataclass <- dataclass(data, model, simS = TRUE)	
	} 
	if(!has.S && K > 1) 
		stop("For starting with sampling parameters 'data' must provide starting allocations.")			
	if(norstud) {
		if(prior@type == "independent") { ## independent prior
			## later regression model ##

			## here only finite mixture ##
			if(length(model@par) == 0) 
				stop("For an independent prior, starting values for the component means have to be provided.")
			has.mu <- "mu" %in% names(model@par)
			if(!has.mu) 
				stop("For an independent prior, starting values for the component means have to be provided.")
		}
		norstudmult <- (model@dist == "normult" || model@dist == "studmult")
		has.logdet <- "logdet" %in% prior@par$sigma  
		if(norstudmult && !has.logdet) {
			has.C <- "C" %in% prior@par$sigma
			if(!has.C) { 
				stop("For an independent prior entry 'C' in 'prior@par$sigma' has to be provided.")
			}
			else { ## if C is there check if array of dimension (r x r x K) 
				if(!is.array(C)) {
					stop("'C' in 'prior@par$sigma' must be an array of dimension (r x r x K).")
				}
				else {
					## check dimensions ##
					dims <- dim(prior@par$sigma$C)
					has.dim <- (dims[1] == r && dims[2] == r && dims[3] == K)
					if(!has.dim) { 
						stop("'C' in 'prior@par$sigma' must be an array of dimension (r x r x K).")
                    }
					logdetC <- array(0, dim = c(r,r,K))
					for(k in 1:K) { ## TODO: check if Cs are given and check priordefine for thi
						logdetC[,,k] <- log(det(prior@par$sigma$C[,,k]))
					}
				}
			}
		} ## end norstudmult
	
	} ## end norstud
	if(model@dist == "cond.poisson") {
		if(!("coef" %in% names(prior@par))) {
			stop("Element 'coef' in slot 'par' of 'prior' is missing. A conditional Poisson mixture needs a coefficient matrix.")
		}
		else {
			coef <- prior@par$coef
			if(!is.matrix(coef) && !is.array(coef))
				stop("Element 'coef' in slot 'par' of 'prior' must be of type 'matrix' or 'array'.")
			if(nrow(coef) != ncol(coef)) 
				stop("ELement 'coef' in slot 'par' of 'prior' must be a quadratic matrix.")
			if(nrow(coef) != K)
				stop("Dimension of element 'coef' in slot 'par' of 'prior' must correspond to the number of components 'K' in 'model'.")
			if(!(all(diag(coef) == 1)))
				stop("Coefficients on the diagnoal of element 'coef' in slot 'par' of 'prior' must be equal to one.")
		}
	} 
			
	######################### MCMC SAMPLING #############################
	if(model@dist == "poisson") {
		## base slots inherited to every derived class ##
		M 		        <- as.integer(mcmc@M)
		ranperm 	    <- mcmc@ranperm
		pars 		    <- list(lambda = array(numeric(), dim = c(M, K)))
		log.mixlik 	    <- array(numeric(), dim = c(M, 1))
		log.mixprior 	<- array(numeric(), dim = c(M, 1))
		## model with fixed indicators ##
		if(model@indicfix || K == 1) { 
			logs 	<- list(mixlik = log.mixlik, mixprior = log.mixprior)
			## model with simple prior ##
			if (!prior@hier) {
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout 	<- new("mcmcoutputfix", M = M, ranperm = ranperm,
								par = pars, log = logs, model = model, 
								prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
				## model output with posterior parameters stored ##
				else {
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.par 	<- list(a = post.a, b = post.b)
					posts  		<- list(par = post.par)
					mcmcout 	<- new("mcmcoutputfixpost", M = M, ranperm = ranperm,
								par = pars, log = logs, post = posts,
								model = model, prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end no hier
			## model with hierarchical prior ##
			else {
				hypers <- list(b = array(numeric(), dim = c(M, 1)))
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout 	<- new("mcmcoutputfixhier", M = M, ranperm = ranperm,
								par = pars, log = logs, hyper = hypers,
								model = model, prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
			
				}
				## model output with posterior parameters stored ##
				else {
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.par 	<- list(a = post.a, b = post.b)
					posts  		<- list(par = post.par)
					mcmcout 	<- new("mcmcoutputfixhierpost", M = M, ranperm = ranperm,
								par = pars, log = logs, hyper = hypers, post = posts,
								model = model, prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end hier 
		} ## end indicfix ##
		## model with simulated indicators ##
		else if(!model@indicfix && K > 1) {			
			log.cdpost 	<- array(numeric(), dim = c(M, 1))
			logs 		<- list(mixlik = log.mixlik, mixprior = log.mixprior, cdpost = log.cdpost)
			weights 	<- array(numeric(), dim = c(M, K))
			entropies 	<- array(numeric(), dim = c(M, 1))
			STm 		<- array(integer(), dim = c(M, 1))
   			Sm 		<- array(integer(), dim = c(N, mcmc@storeS))
			NKm		<- array(integer(), dim = c(M, K))
			clustm 		<- array(integer(), dim = c(N, 1))
			if(mcmc@startpar) {
				Sm[,1] <- as.integer(dataclass$S)
			}
			## model with simple prior ##
			if(!prior@hier) {
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout		<- new("mcmcoutputbase", M = M, ranperm = ranperm,
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, model = model,
								prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
				## model output with posterior parameters stored ##
				else {
					post.weight 	<- array(numeric(), dim = c(M, K))
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.par 	<- list(a = post.a, b = post.b)
					posts  		<- list(par = post.par, weight = post.weight)
					mcmcout 	<- new("mcmcoutputpost", M = M, ranperm = ranperm,
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, post = posts,
								model = model, prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end no hier 
			## model with hierarchical prior ## 
			else {
				hypers 	<- list(b = array(numeric(), dim = c(M, 1)))			
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout	 	<- new("mcmcoutputhier", M = M, ranperm = ranperm, 
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers,
								model = model, prior = prior)
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			## model output with posterior parameters stored ##
				else {
					post.weight 	<- array(numeric(), dim = c(M, K))
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.par 	<- list(a = post.a, b = post.b)
					posts  		<- list(par = post.par, weight = post.weight)
					mcmcout 	<- new("mcmcoutputhierpost", M = M, ranperm = ranperm,
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers, post = posts,
								model = model, prior = prior)		
					.Call("mcmc_poisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end hier
		} ## end no indicfix		
	} ## end model poisson
	
	if(model@dist == "cond.poisson") {
		## base slots inherited to every derived class ##
		M 		<- as.integer(mcmc@M)
		ranperm 	<- mcmc@ranperm
		pars 		<- list(lambda = array(numeric(), dim = c(M, K)))
		log.mixlik 	<- array(numeric(), dim = c(M, 1))
		log.mixprior 	<- array(numeric(), dim = c(M, 1))
		## model with fixed indicators ##
		if(model@indicfix || K == 1) {
			logs 	<- list(mixlik = log.mixlik, mixprior = log.mixprior)
			## model with simple prior ##
			if (!prior@hier) {
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout 	<- new("mcmcoutputfix", M = M, ranperm = ranperm,
								par = pars, log = logs, model = model, 
								prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
				## model output with posterior parameters stored ##
				else {
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.cond	<- array(numeric(), dim	= c(M, K))
					post.par 	<- list(a = post.a, b = post.b, cond = post.cond)
					posts  		<- list(par = post.par)
					mcmcout 	<- new("mcmcoutputfixpost", M = M, ranperm = ranperm,
								par = pars, log = logs, post = posts,
								model = model, prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end no hier
			## model with hierarchical prior ##
			else {
				hypers <- list(b = array(numeric(), dim = c(M, 1)))
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout 	<- new("mcmcoutputfixhier", M = M, ranperm = ranperm,
								par = pars, log = logs, hyper = hypers,
								model = model, prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
			
				}
				## model output with posterior parameters stored ##
				else {
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.cond 	<- array(numeric(), dim	= c(M, K))
					post.par 	<- list(a = post.a, b = post.b, cond = post.cond)
					posts  		<- list(par = post.par)
					mcmcout 	<- new("mcmcoutputfixhierpost", M = M, ranperm = ranperm,
								par = pars, log = logs, hyper = hypers, post = posts,
								model = model, prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end hier 
		} ## end indicfix ##
		## model with simulated indicators ##
		else if(!model@indicfix) {			
			log.cdpost 	<- array(numeric(), dim = c(M, 1))
			logs 		<- list(mixlik = log.mixlik, mixprior = log.mixprior, cdpost = log.cdpost)
			weights 	<- array(numeric(), dim = c(M, K))
			entropies 	<- array(numeric(), dim = c(M, 1))
			STm 		<- array(integer(), dim = c(M, 1))
   			Sm 		<- array(integer(), dim = c(N, mcmc@storeS))
			NKm		<- array(integer(), dim = c(M, K))
			clustm 		<- array(integer(), dim = c(N, 1))
			if(mcmc@startpar) {
				Sm[,1] <- as.integer(dataclass$S)
			}
			## model with simple prior ##
			if(!prior@hier) {
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout		<- new("mcmcoutputbase", M = M, ranperm = ranperm,
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, model = model,
								prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
				## model output with posterior parameters stored ##
				else {
					post.weight 	<- array(numeric(), dim = c(M, K))
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.cond 	<- array(numeric(), dim = c(M, K))
					post.par 	<- list(a = post.a, b = post.b, cond = post.cond)
					posts  		<- list(par = post.par, weight = post.weight)
					mcmcout 	<- new("mcmcoutputpost", M = M, ranperm = ranperm,
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, post = posts,
								model = model, prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end no hier 
			## model with hierarchical prior ## 
			else {
				hypers 	<- list(b = array(numeric(), dim = c(M, 1)))			
				## model output with NO posterior parameters stored ##
				if(!mcmc@storepost) {
					mcmcout	 	<- new("mcmcoutputhier", M = M, ranperm = ranperm, 
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers,
								model = model, prior = prior)
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			## model output with posterior parameters stored ##
				else {
					post.weight 	<- array(numeric(), dim = c(M, K))
					post.a 		<- array(numeric(), dim = c(M, K))
					post.b 		<- array(numeric(), dim = c(M, K))
					post.cond	<- array(numeric(), dim = c(M, K))
					post.par 	<- list(a = post.a, b = post.b, cond = post.cond)
					posts  		<- list(par = post.par, weight = post.weight)
					mcmcout 	<- new("mcmcoutputhierpost", M = M, ranperm = ranperm,
								par = pars, log = logs, weight = weights, entropy = entropies,
								ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers, post = posts,
								model = model, prior = prior)		
					.Call("mcmc_condpoisson_cc", data, model, prior, mcmc, mcmcout, PACKAGE = "finmix")
					return(mcmcout)
				}
			} ## end hier
		} ## end no indicfix		
	} ## end model cond.poisson
} ## end mixturemcmc

