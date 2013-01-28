"mcmcstart" <- function(data, model, varargin) {
	
	## check if mcmc object was given in arguments ##	
	if(nargs() == 2) {
		mcmc <- mcmc()
	}		
	else {
		mcmc <- varargin
	}
	dist <- model@dist
	K <- model@K
	
	## check if data object contains data ##
	has.data <- !all(is.na(data@y))
	has.exposures <- !all(is.na(data@exp))
	if(has.data) {
		if(data@bycolumn) {
			datam <- data@y
			r <- data@r
			if(has.exposures) {
				exp <- data@exp
			}
		}
		else {
			datam <- t(data@y)
			r <- data@r
			if(has.exposures) {
				exp <- t(data@exp)
			}	
		}
	}
	else {
		return("[Error] 'data' object has no data.")
	}

	has.par <- !all(is.na(model@par))
	## check if model object for student-t distributions has a parameter 'df' ##
	if(dist == "student" || dist == "studmult") {
		if(has.par) {
			df.in.model <- "df" %in% names(model@par)
			if(!df.in.model) {	
				new.par <- list(model@par, df = array(10, dim = c(1, K)))
				model@par <- new.par
				validObject(model)	
			}
		}
		else {
			model@par <- list(df = array(10, dim = c(1, K)))
		}
	}

	## check if weights have been already initialized in case of mcmc@startpar == TRUE ##	
	if(K > 1) {
		if(model@indicmod) {
			if(mcmc@starpar && is.na(model@weight)) {
				model@weight <- matrix(1/K, nrow = 1, ncol = K)
			}
		}
	}
	else { ## Markov model, implemented later

	}
		
	## starting values for finite mixtures ## 
	
	## poisson mixtures ##
	if(dist == "poisson") {
		if(!has.par) {
			if(has.exposures) {
				if(K == 1) { ## WARNING: exposures are N x 1!! Change here something in the code
					pm <- max(mean(datam/exp, na.rm = TRUE), 0.1)
					pm <- array(pm, dim = c(1, K))
				}
				else { ## K > 1
					pm <- (mean(datam/exp, na.rm = TRUE)) * exp(runif(K))
					pm <- apply(pm, 2, max, 0.1)
					pm <- array(pm, dim = c(1, K)) 
				}
				model@par <- list(lambda = pm)
			}
			else { ## no exposures (exposures = 1)
				if(K == 1) {
					pm <- max(mean(datam, na.rm = TRUE), 0.1)
					pm <- array(pm, dim = c(1, K))
				}
				else { ## K > 1
					pm <- mean(datam, na.rm) * exp(runif(K))
					pm <- apply(pm, 2, max, 0.1)
					pm <- array(pm, dim = c(1, K))
				}
				model@par <- list(lambda = pm)	
			}
		} 
	}

	## exponential mixtures ##
	else if(dist == "exponential") {
		if(K == 1) {
			pm <- 1/mean(datam, na.rm = TRUE)
			pm <- array(pm, dim = c(1, K))
		}	
		else { ## K > 1
			pm <- exp(runif(K))/mean(datam, na.rm = TRUE)
		}
		model@par <- list(lambda = pm)
	}
	
	## binomial mixtures ##
	else if(dist == "binomial") {
		if(K == 1) {
			pm <- mean(datam, na.rm = TRUE)/model@par$n
			pm <- min(max(pm, 0.1),0.9)
			pm <- array(pm, dim = c(1, K))
		}
		else { ## K > 1
			pm <- mean(datam, na.rm = TRUE)/model@par$n * exp(runif(K))
			pm <- min(max(pm, 0.1), 0.9)
			pm <- array(pm, dim = c(1, K))
		}
		model@par <- list(p = pm, n = model@par$n)
	}

	## univariate normal or student mixtures ##
	else if(dist == "normal" || dist == "student") {
		start.mu <- FALSE
		start.sigma <- FALSE
		
		if(!has.par) {
			start.mu <- TRUE
			start.sigma <- TRUE
		}
		else { ## has already parameters 
			has.mu <- "mu" %in% names(model@par)
			has.sigma <- "sigma" %in% names(model@par)
			if(!has.mu) {
				start.mu <- TRUE
			}
			if(!has.sigma) {
				start.sigma <- TRUE
			}
		}
		
		if(start.mu) {
			if(K == 1) {
				pm <- mean(datam, na.rm = TRUE) 
				pm <- array(pm, dim = c(1, K))
			}
			else { ## K > 1
				pm <- mean(datam, na.rm = TRUE) + sd(datam, na.rm = TRUE) * runif(K)
				pm <- array(pm, dim = c(1, K))
			}
			if(start.sigma) {
				model@par <- list(mu = pm)
			}
			else {
				model@par <- list(mu = pm, model@par)
			}
		}
		
		if(start.sigma) {
		
			pm <- sd(datam, na.rm = TRUE)
			pm <- array(pm, dim = c(1, K))
			model@par <- list(model@par, sigma = pm)
		}
	}

	## multivariate normal mixtures ##
	else if(dist == "normult") {
		start.mu <- FALSE
		start.sigma <- FALSE
		
		if(!has.par) {
			start.mu <- TRUE
			start.sigma <- TRUE
		}
		else {
			has.mu <- "mu" %in% names(model@par)
			has.sigma <- "sigma" %in% names(model@par)
			if(!has.mu) {
				start.mu <- TRUE
			}
			if(!has.sigma) {
				start.sigma <- TRUE
			}
		}
		
		cov.m <- cov(datam)

		if(start.mu) {
			if(K == 1) {
				pm.mu <- apply(datam, 2, mean, na.rm = TRUE)
				pm.mu <- array(t(pm.mu), dim = c(1, K))
			}
			else { ## K > 1
				mean <- apply(datam, 2, mean, na.rm = TRUE)
				pm.mu <- matrix(0, nrow = r, ncol = K)
				for(i in 1:K) {
					pm.mu[,i] = matrix(mean) + chol(cov.m) %*% matrix(runif(K))
				}
			}
			if(start.sigma) {
				model@par <- list(mu = pm.mu)
			}
			else {
				model@par <- list(mu = pm.mu, model@par)
			}
		}
		if(start.sigma) {
			pm.sigma <- array(cov.m, dim = c(r, r, K))
			model@par <- list(model@par, sigma = pm.sigma)
		}
	}	


	## start with classifications ##
	has.S <- !all(is.na(data@S))

	## poisson mixtures ##
	if(dist == "poisson" || dist == "exponential") {
		if(!has.S && K > 1) { ## could possibly produce empty components 
			data@S <- kmeans(datam^.5, centers = K, nstart = K)$cluster
		}
	}

	## binomial mixtures ##
	else if(dist == "binomial") {
		if(!has.S && K > 1) {
			if((max(datam) - min(datam)) > 2 * K) {
				## use k-means to determine a starting classification
				data@S <- kmeans(datam^.5, cemters = K, nstart = K)$cluster
			}
			else {
				## random classification
				data@S <- sample(c(1:K), nrow(datam))
			}
		}
	}

	## univariate normal or student mixtures ##
	else if(dist == "normal" || dist == "student") {
		start.mu <- TRUE
	
		## check if model has parameter mu specified ##
		if(has.par) {
			has.mu <- "mu" %in% names(model@par)
			if(has.mu) {
				start.mu <- FALSE
			}
		}
	
		if(start.mu && K == 1) {
			pm.mu <- mean(datam, na.rm = TRUE)
			pm.mu <- array(na.rm, dim = c(1, K))
			model@par <- list(mu = pm.mu)
		}
	
		else if(has.S && K > 1 && start.mu) {
			if(data@bycolumn) {
				classm <- data@S
			}
			else {
				classm <- t(data@S)
			}
			pm.mu <- array(0, dim = c(1, K))
			for(i in 1:K) {
				pm.mu[i] <- mean(datam[classm == i], na.rm = TRUE) 
			} 
			
		}
	
		else if(!has.S && K > 1) {
			k.means <- kmeans(datam^.5, centers = K, nstart = K)$cluster
			data@S <- k.means$cluster
			if(start.mu) {
				pm.mu <- k.means$centers
				model@par <- list(mu = array(pm.mu, dim = c(1, K)))
			}
		}
	}
	
	## multivariate normal or student-t mixtures ##
	else if(dist == "normult" || dist == "studmult") {
		start.mu <- TRUE

		if(has.par) {
			has.mu <- "mu" %in% names(model@par)
			if(!has.mu) {
				start.mu <- FALSE
			}
		}	
		if(start.mu && K == 1) {
			pm.mu <- apply(datam, 2, mean, na.rm = TRUE)
			pm.mu <- array(na.rm, dim = c(r, K))
			model@par <- list(mu = pm.mu)
		}
	
		else if(has.S && K > 1 && start.mu) {
			if(data@bycolumn) {
				classm <- data@S
			}
			else {
				classm <- t(data@S)
			}
			pm.mu <- array(0, dim = c(r, K))
			for(i in 1:K) {
				pm.mu[i] <- array(apply(datam[classm == i], 2, mean, na.rm = TRUE), dim = c(3, 1)) 
			} 
			
		}
	
		else if(!has.S && K > 1) {
			k.means <- kmeans(datam^.5, centers = K, nstart = K)$cluster
			data@S <- k.means$cluster
			if(start.mu) {
				pm.mu <- k.means$centers
				model@par <- list(mu = array(t(pm.mu), dim = c(r, K)))
			}
		}
	
	}
}

