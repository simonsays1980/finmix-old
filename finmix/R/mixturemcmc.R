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
# along with finmix. If not, see <http://www.gnu.org/licenses/>.

"mixturemcmc" <- function(fdata, model, prior, mcmc) {
    ## Check arguments 
    mcmc <- .check.args.Mixturemcmc(fdata, model, prior, mcmc, nargs())

    ## Default ordering for MCMC: bycolumn
    setBycolumn(fdata)  <- TRUE
    ######################### MCMC SAMPLING #############################
    ## Set the indicators as a default to one for K == 1
    if (model@K == 1) {
        fdata@S <- matrix(1, nrow = fdata@N, ncol = 1)
    }
	if (model@dist == "poisson") {        
        .do.MCMC.Poisson(fdata, model, prior, mcmc)
	} else if (model@dist == "cond.poisson") {
        .do.MCMC.CondPoisson(fdata, model, prior, mcmc)
    }
} ## end mixturemcmc

### Private functions
### These functions are not exported

### Checking
### Check arguments: 'fdata' must contain valid data in @y and in case of 
### starting with sampling the parameters indicators in @S. Further,
### the data in @y must match with the specified distribution in @dist
### of 'model'.
### If it should started with sampling the indicators, 'model' must
### contain valid starting parameters in @par and @weight.
### The 'prior' object must contain valid parameters for the prior
### distribution. 
### Further, if a fixed indicator model is used, @startpar in 'mcmc'
### must be TRUE and @ranperm must be FALSE.
".check.args.Mixturemcmc" <- function(fdata.obj, model.obj,
                                      prior.obj, mcmc.obj, n.args)
{
    ## Check if all arguments are provided
    if (n.args < 4) {
        stop("All arguments must be provided.")
    }
    ## Check if 'fdata' object is valid
    if (class(fdata.obj) != "fdata") {
        stop(paste("Unkown argument. Argument 1 must be an ",
                   "object of class 'fdata'.", sep = ""))
    }
    hasY(fdata.obj, verbose = TRUE)
    ## Check if 'model' was provided:
    if (class(model.obj) != "model") {
        stop(paste("Unknown argument. Argument 2 must be an ",
                   "object of class 'model'.", sep = ""))
    }
    ## Check if 'prior' was provided:
    if (class(prior.obj) != "prior") {
        stop(paste("Unknown argument. Argument 3 must be an ",
                   "object of class 'prior'.", sep = ""))
    }
    ## Check if 'mcmc' was provided:
    if (class(mcmc.obj) != "mcmc") {
        stop(paste("Unkown argument. Argument 4 must be an ",
                   "object of class 'mcmc'.", sep = ""))
    }
    ## Check if @startpar in 'mcmc' object and @indicfix in 
    ## 'model' object match.
    ## For fixed indicator models indicators are not sampled.
    if (model.obj@indicfix && !mcmc.obj@startpar) {
        mcmc.obj@startpar   <- TRUE
    }
    ## Check if @K in 'model' object is one. For a model with
    ## only one component indicators are not sampled.
    if (model.obj@K == 1) {
        mcmc.obj@startpar   <- TRUE
    }
    ## If @startpar in 'mcmc.obj' is TRUE, it should be started 
    ## by sampling the parameters. In this case starting 
    ## indicators must be provided in the 'fdata.obj' object.
    ## If @startpar in 'mcmc.obj' is FALSE it should be started
    ## by sampling the indicators. In this case starting 
    ## parameters must be provided in the 'model.obj' object.
    if (model.obj@K > 1) {
        if (mcmc.obj@startpar) {
            if (!hasS(fdata.obj)) {
                stop(paste("For starting with sampling the parameters ",
                           "the 'fdata' object must contain starting ",
                           "indicator values. See ?mcmcstart for ",
                           "generating valid starting values.", sep = ""))
            }
        } else {
            if (!hasPar(model.obj)) {
                stop(paste("For starting with sampling the indicators ",
                           "the 'model' object must contain starting ",
                           "parameter values. See ?mcmcstart for ",
                           "generating valid starting values.", sep = ""))
            }
            if (!hasWeight(model.obj)) {
                stop(paste("For starting with sampling the indicators ",
                           "the 'model' object must contain starting ",
                           "weight values. See ?mcmcstart for ",
                           "generating valid starting values.", sep = ""))
            }
        }
    }
    ## Check if 'fdata' object and 'model' objects match
    ## Call '.check.fdata.model.Mcmcstart()' from 'mcmcstart.R'.
    .check.fdata.model.Mcmcstart(fdata.obj, model.obj)
    ## Check if 'prior' object is valid    
    if (!hasPriorPar(prior.obj, model.obj)) {
        stop(paste("Slot 'par' in 'prior' object is empty. ",
                   "For MCMC sampling the prior needs fully ",
                   "specified parameters. See ?priordefine for ",
                   "generating valid prior parameters.", sep = ""))
    }
    if (!model.obj@indicfix && model.obj@K > 1) {
        if(!hasPriorWeight(prior.obj, model.obj)) {
            stop(paste("Slot 'weight' of 'prior' object is empty. ",
                       "For MCMC sampling the prior needs specified ",
                       "parameters for the prior of the weights. See ",
                       "?priordefine for generating valid prior ",
                       "parameters.", sep = ""))
        }
    }
    ## Check if @indicfix in 'model' object and 
    ## @ranperm in 'mcmc' object match.
    ## For a fixed indicator model random permutation
    ## sampling is senseless.
    if (model.obj@indicfix && mcmc.obj@ranperm) {
        mcmc.obj@ranperm    <- FALSE
    }
    ## For a model with only one component random permutation
    ## is senseless as well.
    if (model.obj@K == 1 && mcmc.obj@ranperm) {
        mcmc.obj@ranperm    <- FALSE
    }
    return(mcmc.obj)
}
       
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
    } else { ## then check in model 
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
    ## Check for identifiability ##
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

### MCMC
### For each model the MCMC output has to be prepared
### MCMC Poisson: Prepares all data containers for MCMC sampling for
### Poisson mixture models regarding the specifications in 'prior.obj'
### 'model.obj' and 'mcmc.obj'.
".do.MCMC.Poisson" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) 
{
    ## Base slots inherited to every derived class
    K               <- model.obj@K
    N               <- fdata.obj@N
    M 		        <- mcmc.obj@M
    ranperm 	    <- mcmc.obj@ranperm
    burnin          <- mcmc.obj@burnin
    ## Set for MCMC default exposures:
    if (!hasExp(fdata.obj)) {
        fdata.obj@exp   <- matrix(1, nrow = N, ncol = 1)
    }
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
    ## Model with fixed indicators
    if (model.obj@indicfix || K == 1) { 
        logs 	<- list(mixlik = log.mixlik, mixprior = log.mixprior)
        ## Model with simple prior
        if (!prior.obj@hier) {
            ## Model output with NO posterior parameters stored
            if (!mcmc.obj@storepost) {	
                mcmcout 	<- .mcmcoutputfix(M = M, burnin = burnin,
                                              ranperm = ranperm,
                                              par = pars, log = logs,
                                              model = model.obj, 
                                              prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {
            ## Model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputfixpost(M = M, burnin = burnin,
                                                  ranperm = ranperm,
                                                  par = pars, log = logs, 
                                                  post = posts,
                                                  model = model.obj, 
                                                  prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
            ## end no hier
        } else {
        ## Model with hierarchical prior ##
            hypers <- list(b = array(numeric(), dim = c(M, 1)))
            ## Model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout 	<- .mcmcoutputfixhier(M = M,  burnin = burnin,
                                                  ranperm = ranperm,
                                                  par = pars, log = logs, 
                                                  hyper = hypers,
                                                  model = model.obj, 
                                                  prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)			
            } else {
            ## Model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputfixhierpost(M = M, burnin = burnin,
                                                      ranperm = ranperm,
                                                      par = pars, log = logs, 
                                                      hyper = hypers, post = posts,
                                                      model = model.obj, 
                                                      prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
        ## end hier
        }
        ## end indicfix
    } else if (!model.obj@indicfix && K > 1) {			
    ## Model with simulated indicators ##
        log.cdpost 	<- array(numeric(), dim = c(M, 1))
        logs 		<- list(mixlik = log.mixlik, 
                             mixprior = log.mixprior, 
                             cdpost = log.cdpost)
        weights 	<- array(numeric(), dim = c(M, K))
        entropies 	<- array(numeric(), dim = c(M, 1))
        STm 		<- array(integer(), dim = c(M, 1))
        Sm 		    <- array(integer(), dim = c(N, mcmc.obj@storeS))
        NKm		    <- array(integer(), dim = c(M, K))
        clustm 		<- array(integer(), dim = c(N, 1))
        if (!mcmc.obj@startpar) {
            ## First sample for the indicators 
            datac  <- dataclass(fdata.obj, model.obj, simS = TRUE)
            Sm[,1] <- as.integer(datac$S)
        }
        ## Model with simple prior ##
        if (!prior.obj@hier) {
            ## Model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout		<- .mcmcoutputbase(M = M, burnin = burnin,
                                                  ranperm = ranperm,
                                                  par = pars, log = logs, 
                                                  weight = weights, 
                                                  entropy = entropies,
                                                  ST = STm, S = Sm, NK = NKm, 
                                                  clust = clustm, model = model.obj,
                                                  prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            } else {
            ## Model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputpost(M = M, burnin = burnin,
                                               ranperm = ranperm,
                                               par = pars, log = logs, 
                                               weight = weights, 
                                               entropy = entropies,
                                               ST = STm, S = Sm, NK = NKm, 
                                               clust = clustm, post = posts, 
                                               model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            }
        ## end no hier
        } else {
        ## model with hierarchical prior ## 
            hypers 	<- list(b = array(numeric(), dim = c(M, 1)))			
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout	 	<- .mcmcoutputhier(M = M, burnin = burnin,
                                                   ranperm = ranperm, 
                                                   par = pars, log = logs, 
                                                   weight = weights, 
                                                   entropy = entropies,
                                                   ST = STm, S = Sm, NK = NKm, 
                                                   clust = clustm, hyper = hypers,
                                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
               return(mcmcout)
            } else {	
            ## model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputhierpost(M = M, burnin = burnin,
                                                   ranperm = ranperm,
                                                   par = pars, log = logs, 
                                                   weight = weights, 
                                                   entropy = entropies,
                                                   ST = STm, S = Sm, NK = NKm, 
                                                   clust = clustm, hyper = hypers, 
                                                   post = posts,
                                                   model = model.obj, prior = prior.obj)		
                .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            }
        } ## end hier
    } ## end no indicfix		
}

".do.MCMC.CondPoisson" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) 
{
    ## base slots inherited to every derived class ##
    K               <- model.obj@K
    N               <- fdata.obj@N
    M 		        <- mcmc.obj@M
    burnin          <- mcmc.obj@burnin
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
                mcmcout 	<- .mcmcoutputfix(M = M, burnin = burnin,
                                              ranperm = ranperm,
                                              par = pars, log = logs, 
                                              model = model.obj, 
                                              prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputfixpost(M = M, burnin = burnin,
                                                  ranperm = ranperm,
                                                  par = pars, log = logs, 
                                                  post = posts,
                                                  model = model.obj, 
                                                  prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)
            }
            ## end no hier
        } else {
        ## model with hierarchical prior ##
            hypers <- list(b = array(numeric(), dim = c(M, 1)))
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout 	<- .mcmcoutputfixhier(M = M, burnin = burnin,
                                                  ranperm = ranperm,
                                                  par = pars, log = logs, 
                                                  hyper = hypers,
                                                  model = model.obj, 
                                                  prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                return(mcmcout)			
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputfixhierpost(M = M, burnin = burnin,
                                                      ranperm = ranperm,
                                                      par = pars, log = logs, 
                                                      hyper = hypers, 
                                                      post = posts,
                                                      model = model.obj, 
                                                      prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
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
        Sm 		    <- array(integer(), dim = c(N, mcmc.obj@storeS))
        NKm		    <- array(integer(), dim = c(M, K))
        clustm 		<- array(integer(), dim = c(N, 1))
        if (!mcmc.obj@startpar) {
            ## First sample for the indicators 
            datac  <- dataclass(fdata.obj, model.obj, simS = TRUE)
            Sm[,1] <- as.integer(datac$S)
        }
        ## model with simple prior ##
        if (!prior.obj@hier) {
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout		<- .mcmcoutputbase(M = M, burnin = burnin,
                                                  ranperm = ranperm,
                                                  par = pars, log = logs, 
                                                  weight = weights, 
                                                  entropy = entropies,
                                                  ST = STm, S = Sm, NK = NKm, 
                                                  clust = clustm, 
                                                  model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            } else {
            ## model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputpost(M = M, burnin = burnin,
                                               ranperm = ranperm,
                                               par = pars, log = logs, 
                                               weight = weights, 
                                               entropy = entropies,
                                               ST = STm, S = Sm, NK = NKm, 
                                               clust = clustm, post = posts,
                                               model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            }
        ## end no hier
        } else {
        ## model with hierarchical prior ## 
            hypers 	<- list(b = array(numeric(), dim = c(M, 1)))			
            ## model output with NO posterior parameters stored ##
            if (!mcmc.obj@storepost) {
                mcmcout	 	<- .mcmcoutputhier(M = M, burnin = burnin,
                                                   ranperm = ranperm, 
                                                   par = pars, log = logs, 
                                                   weight = weights, 
                                                   entropy = entropies,
                                                   ST = STm, S = Sm, NK = NKm, 
                                                   clust = clustm, hyper = hypers,
                                                   model = model.obj, prior = prior.obj)
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            } else {	
            ## model output with posterior parameters stored ##
                mcmcout 	<- .mcmcoutputhierpost(M = M, burnin = burnin, 
                                                   ranperm = ranperm,
                                                   par = pars, log = logs, 
                                                   weight = weights, 
                                                   entropy = entropies,
                                                   ST = STm, S = Sm, NK = NKm, 
                                                   clust = clustm, hyper = hypers, 
                                                   post = posts,
                                                   model = model.obj, prior = prior.obj)		
                .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj, 
                      mcmc.obj, mcmcout, PACKAGE = "finmix")
                if (mcmc.obj@storeS == 0) {
                    mcmcout@S <- as.array(NA)
                }
                return(mcmcout)
            }
        } ## end hier
    } ## end no indicfix		
}

