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
	if (nargs() < 4) {       
		stop("All arguments must be provided.")
    }

	## check input objects ##
	validObject(data)
    validObject(model)
    validObject(prior)
    validObject(mcmc)

    ## Set global variables ##
	K <- model@K

    ## For fixed indicators random permutation ##
    ## is senseless ##
	if (model@indicfix) {
		ranperm <- FALSE
	}
## --------------------------------------------------------------------------
    ## Check the data ##
	has.data <- !all(is.na(data@y))
	has.S <- !all(is.na(data@S))
	if (has.data) {
		if (data@bycolumn) {
			datam <- data@y
			if(has.S && K == 1) {
				classm <- matrix(1, nrow = nrow(datam), ncol = 1)
			} else if(has.S && K > 1) {
				classm <- data@S
			}
		} else { ## data stored by row
			datam <- t(data@y)
			if (has.S && K == 1) {
				classm <- matrix(1, nrow = nrow(datam), ncol = 1)
			} else if (has.S && K > 1) {
				classm <- t(data@S)
			}
		}
		r <- ncol(datam)
		N <- nrow(datam)
	} else { ## data has no observations
        stop("Observations 'y' of 'data' object 
             are obligatory for MCMC sampling.")
	}
## --------------------------------------------------------------------------   
    ## Check for starting with sampling the parameters ##
    ## or with sampling the allocations ##
	if (mcmc@startpar && !model@indicfix && K > 1) { 
        ## i.e. it should be started by sampling allocations
		if (length(model@par) == 0) {
			stop("For starting with sampling allocations 'model' 
                 must provide starting parameters.")                 
        } ## TODO: else check for each distribution starting parameters.
		if (any(is.na(model@weight)) && model@indicmod == "multinomial") {
			stop("For starting with sampling allocations 'model' 
                 must provide starting weights.")
        }
        ## If all OK: sample allocations.
		dataclass <- dataclass(data, model, simS = TRUE)	
	} else if (!has.S && K > 1) { 
		stop("For starting with sampling the parameters 'data' object 
             must provide starting allocations in 'S'.")			
    }
## --------------------------------------------------------------------------
    ## Check mixture model specific preconditions for ##
    ## MCMC sampling ## 
	if (model@dist == "poisson") {
        data@exp <- .valid.Poisson(data) 
    } else if (model@dist == "cond.poisson") {
        data@exp <- .valid.CondPoisson(data, model, prior)
    } else if (model@dist == "binomial") {
        validObject(model)
        T <- .valid.Reps.Binomial(data, model)
    }
	norstud <- (model@dist == "normal" || model@dist == "normult" 
                || model@dist == "student" || model@dist == "studmult")
    if (norstud) {
        .valid.Norstud(prior, model)
	} ## end norstud
## --------------------------------------------------------------------------
    ######################### MCMC SAMPLING #############################
	if (model@dist == "poisson") {
        .do.MCMC.Poisson(data, model, prior, mcmc, dataclass)
	} else if (model@dist == "cond.poisson") {
        .do.MCMC.CondPoisson(data, model, prior, mcmc, dataclass)
    }
} ## end mixturemcmc

### Private functions
### These functions are not exported

### Validity
### For a Binomial model either the 'data' object
### or the 'model' object must have specified 
### repetitions 'T'. This can be either a 'matrix' 
### object of dimension N x 1 or 1 x 1 (if all 
### repetitions are the same)
".valid.Reps.Binomial" <- function(data, model)
{
    has.reps <- !all(is.na(data@T))
    if (has.reps) {
        if (data@bycolumn) {
            if (nrow(data@T) != N && nrow(data@T) != 1)  {
                stop("Number of repetitions 'T' in 'data' object does 
                     not match number of observations 'N' in 'y'.")
            } else if (nrow(data@T) == N) {
                T <- data@T
            } else { ## dimension of T is 1 x 1
                T <- matrix(data@T[1, 1], nrow = N, ncol = 1)
            }
        } else { ## data stored by row 
            if (ncol(data@T) != N && ncol(data@T) != 1) {
                stop("Number of repetitions 'T' in 'data' object does 
                     not match number of observations 'N' in 'y'.")
            } else if(ncol(data@T) == N) {
                T <- t(data@T)
            } else { ## dimension of T is 1 x 1
                T <- matrix(data@T[1, 1], nrow = N, ncol = 1)	
            }
        }
    } else { ## then check in model ##
        has.reps <- !all(is.na(model@T))
        if (has.reps) {
            if(nrow(model@T) != N && nrow(model@T) != 1) {
                stop("Neither 'data' nor 'model' has correctly 
                     specified repetitions 'T' for binomial model.")
            } else if(nrow(model@T) == N) {
                T <- model@T
            } else { ## dimension of T is 1 x 1 
                T <- matrix(model@T[1, 1], nrow = N, ncol = 1)
            }
        } else {
            stop("Neither 'data' object nor 'model' object has 
                 repetitions 'T' for binomial model specified.")
        }
    }
    ## check for identifiability ##
    ## Reference: Teicher (1961) ##
    rep.occ <- table(T)
    if (dim(unique(T))[1] == 1) {
        if (T[1, 1] < 2 * model@K - 1) {
            warning("This binomial mixture model is not identifiable. 
                    For equal repetitions 'T' it must hold T >= 2K - 1. 
                    See Teicher (1961) for reference.")
        }
    } else {
        if (length(rep.occ) != nrow(T)) {
            if (all(dimnames(rep.occ)$T < rep.occ - 1)) {
                warning("This binomial mixture model is not identifiable.
                        For varying repetitions 'T_i' it must hold T_h
                        > r_h - 1, for unique repetitions 'T_h' and their
                        respective occurences 'r_h'. See Teicher (1961) 
                        for reference")
            } else {
                diff         <- diff(sort(unique(T))) 
                if (any(diff < rep.occ[1:(length(diff))])) {
                    warning("This binomial mixture model is not identifiable. 
                            For varying repetitions 'T_i' it must hold T_h 
                            - T_(h+1) >= r_h for unique repetitions 'T_h' and 
                            respective occurrences 'r_h'. See Teicher (1961) 
                            for reference.")
                }
            }
        }
    }
}

### A conditional Poisson mixture must provide a coefficient 
### matrix. This matrix must be a lower triangular matrix with 
### zeros on its diagonal. Further, the matrix must be quadratic.
".valid.CondPoisson" <- function(data.obj, model.obj, prior.obj) 
{
    N   <- data.obj@N
    K   <- model.obj@K
    if (!("coef" %in% names(prior.obj@par))) {
        stop("Element 'coef' in slot 'par' of 'prior' is missing. A conditional 
             Poisson mixture needs a coefficient matrix.")
    } else {
        coef <- prior.obj@par$coef
        if (!is.matrix(coef) && !is.array(coef)) {
            stop("Element 'coef' in slot 'par' of 'prior' must be 
                 of type 'matrix' or 'array'.")
        }		
        if (nrow(coef) != ncol(coef)) {
            stop("ELement 'coef' in slot 'par' of 'prior' must be 
                 a quadratic matrix.")
        }
        if (nrow(coef) != K) {
            stop("Dimension of element 'coef' in slot 'par' of 'prior' 
                 must correspond to the number of components 'K' in 'model'.")
        }
        if (!(all(diag(coef) == 1))) {
            stop("Coefficients on the diagnoal of element 'coef' in slot 'par' 
                 of 'prior' must be equal to one.")
        }
	}
    has.exposures <- !all(is.na(data.obj@exp))
    if (has.exposures) {
        if (data.obj@bycolumn) {
            if (nrow(data.obj@exp) != N && nrow(data.obj@exp) != 1) {
                stop("Number of exposures 'exp' in 'data' does not 
                     match number of observations in 'y'.")
            }
            else if (nrow(data.obj@exp) == N) {
                exp <- data.obj@exp
            }
            else { ## exp has dimension 1 x 1
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }
        } else { ## data stored by row
            if (ncol(data.obj@exp) != N && ncol(data.obj@exp) != 1) {
                stop("Number of exposures 'exp' in 'data' does not 
                     match number of observations in 'y'.")
            } else if (ncol(data.obj@exp) == N) {
                exp <- t(data.obj@exp)
            } else {
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }
        }
    } else { ## no exposures set default all to 1
        exp <- matrix(1, nrow = N, ncol = 1)
    }   
    return(exp)
}

".valid.Poisson" <- function(data.obj) 
{
    has.exposures   <- !all(is.na(data.obj@exp))
    N               <- data.obj@N
    if (has.exposures) {
        if (data.obj@bycolumn) {
            if (nrow(data.obj@exp) != N && nrow(data.obj@exp) != 1) {
                stop("Number of exposures 'exp' in 'data' does not 
                     match number of observations in 'y'.")
            }
            else if (nrow(data.obj@exp) == N) {
                exp <- data.obj@exp
            }
            else { ## exp has dimension 1 x 1
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }
        } else { ## data stored by row
            if (ncol(data.obj@exp) != N && ncol(data.obj@exp) != 1) {
                stop("Number of exposures 'exp' in 'data' does not 
                     match number of observations in 'y'.")
            } else if (ncol(data.obj@exp) == N) {
                exp <- t(data.obj@exp)
            } else {
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }
        }
    } else { ## no exposures set default all to 1
        exp <- matrix(1, nrow = N, ncol = 1)
    }
    return(exp)
}

".validNorstud" <- function(prior.obj, model.obj) 
{
    if (prior.obj@type == "independent") { ## independent prior
        ## later regression model ##

        ## here only finite mixture ##
        if (length(model.obj@par) == 0) { 
            stop("For an independent prior, starting values for the 
                 parameters have to be provided.")
        }
        has.mu <- "mu" %in% names(model.obj@par)
        if (!has.mu) {
            stop("For an independent prior, starting values for the 
                 component means have to be provided.")
        }
	}
    norstudmult <- (model.obj@dist == "normult" || model.obj@dist == "studmult")
    has.logdet <- "logdet" %in% prior.obj@par$sigma  
    if (norstudmult && !has.logdet) {
        has.C <- "C" %in% prior.obj@par$sigma
        if (!has.C) { 
            stop("For an independent prior entry 'C' in 'prior@par$sigma' has to be provided.")
        } else { ## if C is there check if array of dimension (r x r x K) 
            if (!is.array(C)) {
                stop("'C' in 'prior@par$sigma' must be an array of dimension (r x r x K).")
            } else {
                ## check dimensions ##
                dims <- dim(prior.obj@par$sigma$C)
                has.dim <- (dims[1] == r && dims[2] == r && dims[3] == K)
                if (!has.dim) { 
                    stop("'C' in 'prior@par$sigma' must be an array of dimension (r x r x K).")
                }
                logdetC <- array(0, dim = c(r,r,K))
                for (k in 1:K) { ## TODO: check if Cs are given and check priordefine for thi
                    logdetC[,,k] <- log(det(prior.obj@par$sigma$C[,,k]))
                }
            }
        }
    } ## end norstudmult	
}

### Prepare
### For each model the MCMC output has to be prepared 
".do.MCMC.Poisson" <- function(data.obj, model.obj, prior.obj, mcmc.obj,
                               dataclass) 
{
    ## base slots inherited to every derived class ##
    K               <- model.obj@K
    N               <- data.obj@N
    M 		        <- mcmc.obj@M
    ranperm 	    <- mcmc.obj@ranperm
    pars 		    <- list(lambda = array(numeric(), dim = c(M, K)))
    log.mixlik 	    <- array(numeric(), dim = c(M, 1))
    log.mixprior 	<- array(numeric(), dim = c(M, 1))
    if (mcmc.obj@storepost) {
        post.a 		<- array(numeric(), dim = c(M, K))
        post.b 		<- array(numeric(), dim = c(M, K))
        post.par 	<- list(a = post.a, b = post.b)
        posts  		<- list(par = post.par)
        if (!model.obj@indicfix) {
            posts$weight    <- array(numeric(), dim = c(M, K))
        }
    }
    ## model with fixed indicators ##
    if (model.obj@indicfix || K == 1) { 
        logs 	<- list(mixlik = log.mixlik, mixprior = log.mixprior)
        ## model with simple prior ##
        if (!prior.obj@hier) {
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {	
                mcmcout 	<- new("mcmcoutputfix", M = M, ranperm = ranperm,
                                   par = pars, log = logs, model = model.obj, 
                                   prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputfixpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, post = posts,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
            ## end no hier
        } else {
        ## model with hierarchical prior ##
            hypers <- list(b = array(numeric(), dim = c(M, 1)))
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout 	<- new("mcmcoutputfixhier", M = M, ranperm = ranperm,
                                   par = pars, log = logs, hyper = hypers,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)			
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputfixhierpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, hyper = hypers, post = posts,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        ## end hier
        }
        ## end indicfix
    } else if (!model.obj@indicfix && K > 1) {			
    ## model with simulated indicators ##
        log.cdpost 	<- array(numeric(), dim = c(M, 1))
        logs 		<- list(mixlik = log.mixlik, mixprior = log.mixprior, cdpost = log.cdpost)
        weights 	<- array(numeric(), dim = c(M, K))
        entropies 	<- array(numeric(), dim = c(M, 1))
        STm 		<- array(integer(), dim = c(M, 1))
        Sm 		<- array(integer(), dim = c(N, mcmc.obj@storeS))
        NKm		<- array(integer(), dim = c(M, K))
        clustm 		<- array(integer(), dim = c(N, 1))
        if (mcmc.obj@startpar) {
            Sm[,1] <- as.integer(dataclass$S)
        }
        ## model with simple prior ##
        if (!prior.obj@hier) {
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout		<- new("mcmcoutputbase", M = M, ranperm = ranperm,
                                      par = pars, log = logs, weight = weights, entropy = entropies,
                                      ST = STm, S = Sm, NK = NKm, clust = clustm, model = model.obj,
                                      prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, weight = weights, entropy = entropies,
                                   ST = STm, S = Sm, NK = NKm, clust = clustm, post = posts,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        ## end no hier
        } else {
        ## model with hierarchical prior ## 
            hypers 	<- list(b = array(numeric(), dim = c(M, 1)))			
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout	 	<- new("mcmcoutputhier", M = M, ranperm = ranperm, 
                                       par = pars, log = logs, weight = weights, entropy = entropies,
                                       ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers,
                                       model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {	
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputhierpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, weight = weights, entropy = entropies,
                                   ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers, post = posts,
                                   model = model.obj, prior = prior.obj)		
                .Call("mcmc_poisson_cc", data.obj, model.obj, prior.obj, mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        } ## end hier
    } ## end no indicfix		
}

".do.MCMC.CondPoisson" <- function(data.obj, model.obj, prior.obj, 
                                   mcmc.obj, dataclass) 
{
    ## base slots inherited to every derived class ##
    K               <- model.obj@K
    N               <- data.obj@N
    M 		        <- mcmc.obj@M
    ranperm 	    <- mcmc.obj@ranperm
    pars 		    <- list(lambda = array(numeric(), dim = c(M, K)))
    log.mixlik 	    <- array(numeric(), dim = c(M, 1))
    log.mixprior 	<- array(numeric(), dim = c(M, 1))
    if (mcmc.obj@storepost) {
        post.a 		<- array(numeric(), dim = c(M, K))
        post.b 		<- array(numeric(), dim = c(M, K))
        post.par 	<- list(a = post.a, b = post.b)
        posts  		<- list(par = post.par)
        if (!model.obj@indicfix) {
            posts$weight    <- array(numeric(), dim = c(M, K))
        }
    }
    ## model with fixed indicators ##
    if (model.obj@indicfix || K == 1) { 
        logs 	<- list(mixlik = log.mixlik, mixprior = log.mixprior)
        ## model with simple prior ##
        if (!prior.obj@hier) {
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {	
                mcmcout 	<- new("mcmcoutputfix", M = M, ranperm = ranperm,
                                   par = pars, log = logs, model = model.obj, 
                                   prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputfixpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, post = posts,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
            ## end no hier
        } else {
        ## model with hierarchical prior ##
            hypers <- list(b = array(numeric(), dim = c(M, 1)))
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout 	<- new("mcmcoutputfixhier", M = M, ranperm = ranperm,
                                   par = pars, log = logs, hyper = hypers,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)			
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputfixhierpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, hyper = hypers, post = posts,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        ## end hier
        }
        ## end indicfix
    } else if (!model.obj@indicfix && K > 1) {			
    ## model with simulated indicators ##
        log.cdpost 	<- array(numeric(), dim = c(M, 1))
        logs 		<- list(mixlik = log.mixlik, mixprior = log.mixprior, cdpost = log.cdpost)
        weights 	<- array(numeric(), dim = c(M, K))
        entropies 	<- array(numeric(), dim = c(M, 1))
        STm 		<- array(integer(), dim = c(M, 1))
        Sm 		<- array(integer(), dim = c(N, mcmc.obj@storeS))
        NKm		<- array(integer(), dim = c(M, K))
        clustm 		<- array(integer(), dim = c(N, 1))
        if (mcmc.obj@startpar) {
            Sm[,1] <- as.integer(dataclass$S)
        }
        ## model with simple prior ##
        if (!prior.obj@hier) {
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout		<- new("mcmcoutputbase", M = M, ranperm = ranperm,
                                      par = pars, log = logs, weight = weights, entropy = entropies,
                                      ST = STm, S = Sm, NK = NKm, clust = clustm, model = model.obj,
                                      prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, weight = weights, entropy = entropies,
                                   ST = STm, S = Sm, NK = NKm, clust = clustm, post = posts,
                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        ## end no hier
        } else {
        ## model with hierarchical prior ## 
            hypers 	<- list(b = array(numeric(), dim = c(M, 1)))			
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout	 	<- new("mcmcoutputhier", M = M, ranperm = ranperm, 
                                       par = pars, log = logs, weight = weights, entropy = entropies,
                                       ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers,
                                       model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {	
            ## model output with posterior parameters stored ##
                mcmcout 	<- new("mcmcoutputhierpost", M = M, ranperm = ranperm,
                                   par = pars, log = logs, weight = weights, entropy = entropies,
                                   ST = STm, S = Sm, NK = NKm, clust = clustm, hyper = hypers, post = posts,
                                   model = model.obj, prior = prior.obj)		
                .Call("mcmc_condpoisson_cc", data.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        } ## end hier
    } ## end no indicfix		
}

