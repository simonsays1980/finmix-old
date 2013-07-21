### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

".setUp.data" <- function(withInd = FALSE) {
        ## Get path ##
        data.path <- paste(getwd(), 
			"/../data/poisson.data.csv", sep= "")
        data <- read.csv(data.path, header = FALSE, sep = ",")
        if(withInd) {
                ind.path <- paste(getwd(), 
				"/../data/poisson.ind.csv", sep = "")
                ind <- read.csv(ind.path, header = FALSE, sep = ",")
                data <- data(y. = as.matrix(data), S. = as.matrix(ind), type. = "discrete",
                                r. = 1, N. = nrow(data), sim. = TRUE,
                                bycolumn. = TRUE)
                return(data)
        }
        else {
                data <- data(y. = as.matrix(data), type. = "discrete", r. = 1,
                                N. = nrow(data), sim. = TRUE,
                                bycolumn. = TRUE)
                return(data)
        }
}

".setUp.model" <- function() {
        model <- model(dist. = "poisson", K. = 2)
        return(model)
}

".setUp.mcmc" <- function() {
        mcmc <- mcmc(burnin. = 0, M. = 100, startpar. = FALSE,
                        storeS. = 2, storepost. = FALSE,
                        ranperm. = FALSE)
        return(mcmc)
}

### --- Test functions --- ###

"test.subseq.mcmcoutputfix.poisson" <- function() {
      ## --- Check K = 1 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 1
      setIndicfix(model)    <- TRUE
      prior                 <- prior(hier = FALSE)
      prior                 <- priordefine(data, model, varargin = prior)
      mcmc                  <- .setUp.mcmc()
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      sub.index             <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, 
                                      ncol = 1) 
      mcmcout.sub           <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout.sub@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5] == mcmcout.sub@par$lambda), "check4")
      ## Test exception ## 
      sub.index     <- matrix(seq(1:50) < 6, nrow = 50, ncol = 1)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
      ## --- Check for K = 2 --- ##
      set.seed(0)
      setK(model)   <- 2
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check8")
      ## --- Check for K = 3 --- ##
      set.seed(0)
      setK(model)   <- 3
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check12")
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputfixhier.poisson" <- function() {
      ## --- Check K = 1 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 1
      setIndicfix(model)    <- TRUE
      prior                 <- priordefine(data, model)
      mcmc                  <- .setUp.mcmc()
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      sub.index             <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, 
                                      ncol = 1) 
      mcmcout.sub           <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout.sub@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout@hyper$b[1:5] == mcmcout.sub@hyper$b))
      ## Test exception ## 
      sub.index     <- matrix(seq(1:50) < 6, nrow = 50, ncol = 1)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
      ## --- Check for K = 2 --- ##
      set.seed(0)
      setK(model)   <- 2
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check7")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout@hyper$b[1:5] == mcmcout.sub@hyper$b))
      ## --- Check for K = 3 --- ##
      set.seed(0)
      setK(model)   <- 3
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check10")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout@hyper$b[1:5] == mcmcout.sub@hyper$b))
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

