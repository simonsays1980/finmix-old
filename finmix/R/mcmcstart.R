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
# along with finmix.  If not, see <http://www.gnu.org/licenses/>.

"mcmcstart" <- function(data, model, varargin) 
{
    K       <- model@K
    dist    <- model@dist
	## Check if mcmc object was given in arguments
	if(nargs() == 2) {
		mcmc <- mcmc()
	}		
	else {
		mcmc <- varargin
	}
    ## If it is started by sampling the allocations it needs start 
    ## parameters. 
    if (mcmc@startpar) {
        ## Check if model object for student-t distributions has 
        ## a parameter 'df'.
        if (dist %in% c("student", "studmult")) {
            .mcmcstart.Student.Df(model.obj)
        }
    	## Check if weights have been already initialized
        if (K > 1) {
    		if (model@indicmod == "multinomial") {
                model <- .mcmcstart.Multinomial.Model(model, mcmc)
	    	} else { ## Markov model, implemented later.
            
            }
        }
        if (dist %in% c("poisson", "cond.poisson")) {
            model <- .mcmcstart.Poisson.Model(data, model, mcmc)
        } else if (dist == "exponential") {
            model <- .mcmcstart.Exponential.Model(data, model, mcmc)
        } else if (dist == "binomial") {
            model <- .mcmcstart.Binomial.Model(data, model, mcmc)
        } else if (dist == "normal" || dist == "student") {
            model <- .mcmcstart.Norstud.Model(data, model, mcmc)
        } else if (dist %in% c("normult", "studmult")) {
            model <- .mcmcstart.Norstudmult.Model(data, model, mcmc)
        }
    } else { ## Start by sampling the parameters (default)
        if (K > 1) {
            if (dist %in% c("poisson", "exponential")) {
                data <- .mcmcstart.Ind.Poisson(data, model)
            } else if (dist == "binomial") {
                data <- .mcmcstart.Ind.Binomial(data, model)
            } else if(dist %in% c("normal", "normult", 
                                  "student", "studmult")) {
                data <- .mcmcstart.Ind.Norstud(data, model)
            }        
        }
    }
	olist <- list(data = data, model = model, mcmc = mcmc)
	return(olist)
}

### Private functions.
### These functions are not exported.
".mcmcstart.Data" <- function(data.obj) 
{
    if(data.obj@bycolumn) {
        datam <- data.obj@y
    } else {
        datam <- t(data@y)
    }
    return(datam)
}

".mcmcstart.Exp" <- function(data.obj) 
{
    r           <- data.obj@r
    N           <- data.obj@N
    has.exp     <- (length(data.obj@exp) > 0) 
    if (has.exp) {
        if (data.obj@bycolumn) {
            if (nrow(data.obj@exp) != N && nrow(data.obj@exp) != 1) {
                stop(paste("Dimension of slot 'exp' of 'data' object",
                           "does not match dimension of slot 'y' of",
                           "'data' object."), 
                     sep = "")
            } else if (nrow(data.obj@exp) == N) {
                exp <- data.obj@exp
            } else {
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }
        } else {
            if (ncol(data.obj@exp) != N && ncol(data.obj@exp) != 1) {
                stop(paste("Dimension of slot 'exp' of 'data' object",
                           "does not match dimension of slot 'y' of",
                           "'data' object."),
                     sep = "")
            } else if (ncol(data.obj@exp) == N) {
                exp <- t(data.obj@exp)
            } else {
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }           
        }
    } else {
        exp <- matrix(1, nrow = N, ncol = 1)
    }
    return(exp)
}

".mcmcstart.Multinomial.Model" <- function(model.obj, mcmc.obj)
{
    if (all(is.na(model.obj@weight))) {
        model@weight <- matrix(1/K, nrow = 1, ncol = K)
    }
    return(model.obj)
}

".mcmcstart.Poisson.Model" <- function(data.obj, model.obj, 
                                       mcmc.obj) 
{
    K           <- model.obj@K
    has.par     <- (length(model.obj@par) > 0)
    has.exp     <- !all(is.na(data.obj@exp))                       
    if (mcmc.obj@startpar && !has.par) {
        if (has.exp) {
            if (K == 1) { 
                pm <- max(mean(data.obj@y/exp, na.rm = TRUE), 0.1)
                pm <- array(pm, dim = c(1, K))
            } else { ## K > 1
                pm <- (mean(data.obj@y/exp, na.rm = TRUE)) * exp(runif(K))
                pm <- pmax(pm, 0.1)
                pm <- array(pm, dim = c(1, K)) 
            }
            model.obj@par <- list(lambda = pm)
        }
        else { ## no exposures (exposures = 1)
            if(K == 1) {
                pm <- max(mean(data.obj@y, na.rm = TRUE), 0.1)
                pm <- array(pm, dim = c(1, K))
            } else { ## K > 1
                pm <- mean(data.obj@y, na.rm = TRUE) * exp(runif(K))
                pm <- pmax(pm, 0.1)
                pm <- array(pm, dim = c(1, K))
            }
            model.obj@par <- list(lambda = pm)	
        }
    }
    return(model.obj)
}

".mcmcstart.Exponential.Model" <- function(data.obj, model.obj,
                                           mcmc.obj)
{
    K               <- model.obj@K
    has.par         <- (length(model.obj@par) > 0)
    if (mcmc.obj@startpar && !has.par) {
        if (K == 1) {
            pm <- 1/mean(data.obj@y, na.rm = TRUE)
        } else { ## K > 1
            pm <- exp(runif(K))/mean(data.obj@y, na.rm = TRUE)
        }
        model.obj@par <- list(lambda = pm)
    }
    return(model.obj)
}

".mcmcstart.Binomial.Model" <- function(data.obj, model.obj,
                                        mcmc.obj)
{
    K           <- model.obj@K
    has.par     <- (length(model.obj@par) > 0)
    if (mcmc.obj@startpart && !has.par) {
        if (K == 1) {
            pm <- mean(data.obj@y, na.rm = TRUE)/model.obj@par$n
            pm <- pmin(pmax(pm, 0.1),0.9)
        } else { ## K > 1
            pm <- mean(data.obj@y, na.rm = TRUE)/data.obj@T * exp(runif(K))
            pm <- pmin(pmax(pm, 0.1), 0.9)
        }
        model.obj@par <- list(p = pm)
    }
    return(model.obj)
}

".mcmcstart.Norstud.Model" <- function(data.obj, model.obj,
                                       mcmc.obj)
{
    K           <- model.obj@K
    has.par     <- (length(model.obj@par) > 0)
    if(mcmc@startpar) {
        start.mu    <- FALSE
        start.sigma <- FALSE
	    if(!has.par) {
				start.mu    <- TRUE
				start.sigma <- TRUE
        } else { ## has already parameters 
            has.mu      <- "mu" %in% names(model.obj@par)
            has.sigma   <- "sigma" %in% names(model.obj@par)
            if(!has.mu) {
                start.mu    <- TRUE
            }
            if(!has.sigma) {
                start.sigma <- TRUE
            }
        }		
        if(start.mu) {
            if(K == 1) {
                pm <- mean(data.obj@y, na.rm = TRUE) 
            } else { ## K > 1
                pm <- mean(data.obj@y, na.rm = TRUE) + 
                        sd(data.obj@y, na.rm = TRUE) * runif(K)
                pm <- matrix(pm, nrow = 1, ncol = K)
            }
            if(start.sigma) {
                model.obj@par <- list(mu = pm)
            } else {
                model.obj@par <- list(mu = pm, model.obj@par)
            }
        }
        if(start.sigma) {			
            pm                  <- sd(data.obj@y, na.rm = TRUE)
            pm                  <- matrix(pm, nrow = 1, ncol = K)
            model.obj@par       <- list(model.obj@par, sigma = pm)
        }
    }
    return(model.obj)
}

".mcmcstart.Norstudmult.Model" <- function(data.obj, model.obj,
                                       mcmc.obj)
{
    K               <- model.obj@K
    r               <- model.obj@r
    has.par         <- (length(model.obj@par) > 0)
    datam           <- .mcmcstart.Data(data.obj)
    ## Check if parameters are already provided ##
    start.mu    <- FALSE
    start.sigma <- FALSE
    if (!has.par) {
        start.mu    <- TRUE
        start.sigma <- TRUE
    } else {
        has.mu      <- "mu" %in% names(model.obj@par)
        has.sigma   <- "sigma" %in% names(model.obj@par)
        if (!has.mu) {
            start.mu    <- TRUE
        }
        if (!has.sigma) {
            start.sigma <- TRUE
        }
    }			
    cov.m <- cov(datam)	
    if (start.mu) {
        if (K == 1) {
            pm.mu   <- apply(datam, 2, mean, na.rm = TRUE)
            pm.mu   <- array(pm.mu, dim = c(1, K))
        }
        else { ## K > 1
            mean    <- apply(datam, 2, mean, na.rm = TRUE)
            pm.mu   <- matrix(0, nrow = r, ncol = K)
            for(i in 1:K) {
                pm.mu[,i] <-  matrix(mean) + t(chol(cov.m)) %*% matrix(runif(K))
            }
        }
        if (!has.par) {
            model.obj@par <- list(mu = pm.mu)
        }
        else {
            model.obj@par$mu <- pm.mu
        }
    }
    if (start.sigma) {
        model.obj@par$sigma <- array(cov.m, dim = c(r, r, K))
    }
    return(model.obj)
}

".mcmcstart.Student.Df" <- function(model.obj)
{
    has.par     <- (length(model.obj@par) > 0)
    if (has.par) {
        has.df  <- "df" %in% names(model.obj@par)
        if (!df.in.model) {	
            model.obj@pari$df <- array(10, dim = c(1, K))
            validObject(model.obj)	
        }			
    } else {
        model@par <- list(df = array(10, dim = c(1, K)))
    }
    return(model.obj)
}

### Functions for generating start values for indicators
".mcmcstart.valid.Ind" <- function(data.obj)
{
    N       <- data.obj@N
    has.S   <- (length(data.obj@S) > 0)
    if (has.S) {
        if (data.obj@bycolumn) {
            if (nrow(data.obj@S) != N || ncol(data.obj@S) != 1) {
                stop(paste("Dimension of slot 'S' in 'data'",
                           "object does not match dimension",
                           "of slot 'y' in 'data' object."), sep ="")
            }
        } else {
            if (ncol(data.obj@S) != N || nrow(data.obj@S) != 1) {
                stop(paste("Dimension of slot 'S' in 'data'",
                           "object does not match dimension",
                           "of slot 'y' in 'data' object."), sep = "")
            }
        }
        return(TRUE)
    } else {
        return(FALSE)
    }
}

".mcmcstart.Ind.Poisson" <- function(data.obj, model.obj)
{
    K       <- model.obj@K
    has.S   <- .mcmcstart.valid.Ind(data.obj)
    datam   <- .mcmcstart.Data(data.obj) 
    if (!has.S) {
        if (data.obj@bycolumn) {
            data.obj@S <- matrix(kmeans(datam^.5, 
                                        centers = K, 
                                        nstart = K)$cluster)
        } else {
            data.obj@S <-t( matrix(kmeans(datam^.5, 
                                        centers = K, 
                                        nstart = K)$cluster))
        }
    }
    return(data.obj)
}

".mcmcstart.Ind.Binomial" <- function(data.obj, model.obj) 
{
    K           <- model.obj
    N           <- data.obj
    has.S       <- .mcmcstart.valid.Ind(data.obj)
    datam       <- .mcmcstart.Data(data.obj)
    if((max(datam) - min(datam)) > 2 * K) {
        ## use k-means to determine a starting classification
        if (data.obj@bycolumn) {
            data.obj@S <- as.matrix(kmeans(datam^.5, 
                                           centers = K, 
                                           nstart = K)$cluster)
        } else {
            data.obj@S <- t(as.matrix(kmeans(datam^.5, 
                                             centers = K,
                                             nstart = K)$cluster))
        }
    } else {
        ## random classification
        if (data.obj@bycolumn) {
            data.obj@S  <- as.matrix(sample(c(1:K), N, 
                                            replace = TRUE))
        } else {
            data.obj@S  <- t(as.matrix(sample(c(1:K), N, 
                                              replace = TRUE)))
        }
    }
    return(data.obj)
}

".mcmcstart.Ind.Norstud" <- function(data.obj, model.obj) 
{
    K           <- model.obj@K
    has.S       <- .mcmcstart.valid.Ind(data.obj)
    datam       <- .mcmcstart.Data(data.obj)
    if (has.S) {
        return(data.obj)
    } else {
        if (data.obj@bycolumn) {
            data.obj@S  <- as.matrix(kmeans(datam^.5, 
                                            centers = K, 
                                            nstart = K)$cluster)
        } else {
            data.obj@S  <- t(as.matrix(kmeans(datam^.5, 
                                            centers = K,
                                            nstart = K)$cluster))
        }
        return(data.obj)
    }
}
