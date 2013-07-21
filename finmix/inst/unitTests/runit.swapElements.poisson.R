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

### --- Test [[Rcpp::export]] functions --- ###

"test.swap_cc" <- function() {
    set.seed(0)
    values          <- matrix(rnorm(20), nrow = 10, ncol = 2) 
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 10, ncol = 2, byrow = TRUE)
    values.perm     <- swap_cc(values, perm.index) 
    ## Test cases ##
    checkEquals(nrow(values), nrow(values.perm))
    checkEquals(ncol(values), ncol(values.perm))
    checkTrue(!any(values == values.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 1, ncol = 2, byrow = TRUE)
    checkException(swap_cc(values, perm.index), silent = TRUE)
    ## --- Check for K = 3 --- ##
    values          <- matrix(rnorm(30), nrow = 10, ncol = 3)
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
    values.perm     <- swap_cc(values, perm.index)
    ## Test cases ##
    checkEquals(nrow(values), nrow(values.perm))
    checkEquals(ncol(values), ncol(values.perm))
    checkTrue(!any(values == values.perm), "check7")
}

"test.swapInd_cc" <- function() {
    set.seed(0) 
    indicator       <- matrix(sample(c(1,2), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 1, ncol = 2, byrow = TRUE)
    indicator.perm  <- swapInd_cc(indicator, perm.index)
    ## Test cases ##
    checkEquals(nrow(indicator.perm), nrow(indicator))
    checkEquals(ncol(indicator.perm), ncol(indicator))
    checkTrue(!any(indicator == indicator.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(c(2,1), nrow = 2, ncol = 2, byrow = TRUE)
    checkException(swapInd_cc(indicator, perm.index), silent = TRUE)
    ## --- Check K = 3 --- ##
    set.seed(0)
    indicator       <- matrix(sample(c(1, 2, 3), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 1, ncol = 3, byrow = TRUE)
    indicator.perm  <- swapInd_cc(indicator, perm.index) 
    ## Test cases ##
    checkEquals(nrow(indicator), nrow(indicator.perm))
    checkEquals(ncol(indicator), ncol(indicator.perm))
    checkTrue(!any(indicator == indicator.perm), "check7")
}

"test.swapST_cc" <- function() {
    set.seed(0)
    indicator       <- matrix(sample(c(1,2), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 10, ncol = 2, byrow = TRUE)
    indicator.perm  <- swapST_cc(indicator, perm.index)
    ## Test cases ##
    checkEquals(nrow(indicator), nrow(indicator.perm))
    checkEquals(ncol(indicator), ncol(indicator.perm))
    checkTrue(!any(indicator == indicator.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 2, ncol = 2, byrow = TRUE)
    checkException(swapST_cc(indicator, perm.index), silent = TRUE)
    ## --- Check for K = 3 --- ##
    set.seed(0) 
    indicator       <- matrix(sample(c(1, 2, 3), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
    indicator.perm  <- swapST_cc(indicator, perm.index)
    ## Test cases ##
    checkEquals(nrow(indicator), nrow(indicator.perm))
    checkEquals(ncol(indicator), ncol(indicator.perm))
    checkTrue(!any(indicator == indicator.perm), "check7")
}


"test.swapElements.mcmcoutputfix.poisson" <- function() {
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
      perm.index            <- matrix(as.integer(2), nrow = 100, ncol = 1) 
      mcmcout.perm          <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout.perm@par$lambda), nrow(mcmcout@par$lambda))
      checkEquals(ncol(mcmcout.perm@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      ## Test exception ## 
      perm.index    <- matrix(as.integer(2), nrow = 100, ncol = 2)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
      ## --- Check for K = 2 --- ##
      setK(model)   <- 2
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check7")
      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check10")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)


}

"test.swapElements.mcmcoutputfixpost.poisson" <- function() {
      ## --- Check K = 1 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 1
      setIndicfix(model)    <- TRUE
      prior                 <- prior(hier = FALSE)
      prior                 <- priordefine(data, model, varargin = prior)
      mcmc                  <- .setUp.mcmc()
      setStorepost(mcmc)    <- TRUE
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      perm.index            <- matrix(as.integer(2), nrow = 100, ncol = 1)
      mcmcout.perm          <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout.perm@par$lambda), nrow(mcmcout@par$lambda))
      checkEquals(ncol(mcmcout.perm@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout@post$par$b == mcmcout.perm@post$par$b), "check8")
      checkTrue(all(mcmcout@post$par$a == mcmcout.perm@post$par$a), "check9")
      ## Test exception ## 
      perm.index    <- matrix(as.integer(2), nrow = 100, ncol = 2)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
      ## --- Check for K = 2 --- ##
      set.seed(0)
      setK(model)   <- 2
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check12")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(!any(mcmcout@post$par$b == mcmcout.perm@post$par$b), "check17")
      checkTrue(!any(mcmcout@post$par$a == mcmcout.perm@post$par$a), "check18")
      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check21")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(!any(mcmcout@post$par$b == mcmcout.perm@post$par$b), "check17")
      checkTrue(!any(mcmcout@post$par$a == mcmcout.perm@post$par$a), "check18")
     ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}

"test.swapElements.mcmcoutputfixhier.poisson" <- function() {
      ## --- Check K = 1 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 1
      setIndicfix(model)    <- TRUE
      prior                 <- priordefine(data, model)
      mcmc                  <- .setUp.mcmc()
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      perm.index            <- matrix(as.integer(2), nrow = 100, ncol = 1)
      mcmcout.perm          <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout.perm@par$lambda), nrow(mcmcout@par$lambda))
      checkEquals(ncol(mcmcout.perm@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      ## Test exception ## 
      perm.index    <- matrix(as.integer(2), nrow = 100, ncol = 2)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
      ## --- Check for K = 2 --- ##
      set.seed(0)
      setK(model)   <- 2
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check12")
      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check21")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}

"test.swapElements.mcmcoutputfixhierpost.poisson" <- function() {
      ## --- Check K = 1 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 1
      setIndicfix(model)    <- TRUE
      prior                 <- priordefine(data, model)
      mcmc                  <- .setUp.mcmc()
      setStorepost(mcmc)    <- TRUE
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      perm.index            <- matrix(as.integer(2), nrow = 100, ncol = 1)
      mcmcout.perm          <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout.perm@par$lambda), nrow(mcmcout@par$lambda))
      checkEquals(ncol(mcmcout.perm@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout@post$par$b == mcmcout.perm@post$par$b), "check8")
      checkTrue(all(mcmcout@post$par$a == mcmcout.perm@post$par$a), "check9")
      ## Test exception ## 
      perm.index    <- matrix(as.integer(2), nrow = 100, ncol = 2)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
      ## --- Check for K = 2 --- ##
      set.seed(0)
      setK(model)   <- 2
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check12")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(!any(mcmcout@post$par$b == mcmcout.perm@post$par$b), "check17")
      checkTrue(!any(mcmcout@post$par$a == mcmcout.perm@post$par$a), "check18")
      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check21")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(!any(mcmcout@post$par$b == mcmcout.perm@post$par$b), "check17")
      checkTrue(!any(mcmcout@post$par$a == mcmcout.perm@post$par$a), "check18")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}

"test.swapElements.mcmcoutputind.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data          <- .setUp.data(withInd = TRUE)
      model         <- .setUp.model()
      setK(model)   <- 2
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmc          <-.setUp.mcmc()
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check6")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check9")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check12")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check15")

      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check16")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check19")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check22")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check25")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check28")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}



"test.swapElements.mcmcoutputindpost.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)   <- 2
      prior                 <- prior(hier = FALSE)
      prior                 <- priordefine(data, model, varargin = prior)
      mcmc                  <- .setUp.mcmc()
      setStorepost(mcmc)    <- TRUE
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(any(mcmcout@post$par$b != mcmcout.perm@post$par$b), "check9")
      checkTrue(any(mcmcout@post$par$a != mcmcout.perm@post$par$a), "check10")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check13")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check16")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check19")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check22")
      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check25")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(any(mcmcout@post$par$b != mcmcout.perm@post$par$b), "check30")
      checkTrue(any(mcmcout@post$par$a != mcmcout.perm@post$par$a), "check31")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check34")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check37")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check40")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check43")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}

"test.swapElements.mcmcoutputhier.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data          <- .setUp.data(withInd = TRUE)
      model         <- .setUp.model()
      setK(model)   <- 2
      prior         <- priordefine(data, model)
      mcmc          <-.setUp.mcmc()
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check6")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check9")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check12")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check15")

      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check16")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check19")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check22")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check25")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check28")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}



"test.swapElements.mcmcoutputindhierpost.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 2
      prior                 <- priordefine(data, model)
      mcmc                  <- .setUp.mcmc()
      setStorepost(mcmc)    <- TRUE
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2,1)), nrow = 100, ncol = 2, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check3")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(any(mcmcout@post$par$b != mcmcout.perm@post$par$b), "check9")
      checkTrue(any(mcmcout@post$par$a != mcmcout.perm@post$par$a), "check10")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check13")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check16")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check19")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check22")
      ## --- Check for K = 3 --- ##
      setK(model)   <- 3
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 100, ncol = 3, byrow = TRUE)
      mcmcout.perm  <- swapElements(mcmcout, perm.index)
      ## Test cases ##
      checkEquals(nrow(mcmcout@par$lambda), nrow(mcmcout.perm@par$lambda))
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.perm@par$lambda))
      checkTrue(!any(mcmcout@par$lambda == mcmcout.perm@par$lambda), "check25")
      checkEquals(nrow(mcmcout.perm@post$par$a), nrow(mcmcout@post$par$a))
      checkEquals(ncol(mcmcout.perm@post$par$a), ncol(mcmcout@post$par$a))
      checkEquals(nrow(mcmcout.perm@post$par$b), nrow(mcmcout@post$par$b))
      checkEquals(ncol(mcmcout.perm@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(any(mcmcout@post$par$b != mcmcout.perm@post$par$b), "check30")
      checkTrue(any(mcmcout@post$par$a != mcmcout.perm@post$par$a), "check31")
      checkEquals(nrow(mcmcout.perm@weight), nrow(mcmcout@weight))
      checkEquals(ncol(mcmcout.perm@weight), ncol(mcmcout@weight))
      checkTrue(!any(mcmcout@weight == mcmcout.perm@weight), "check34")
      checkEquals(nrow(mcmcout.perm@ST), nrow(mcmcout@ST))
      checkEquals(ncol(mcmcout.perm@ST), ncol(mcmcout@ST))
      checkTrue(!any(mcmcout@ST == mcmcout.perm@ST), "check37")
      checkEquals(nrow(mcmcout.perm@S), nrow(mcmcout@S))
      checkEquals(ncol(mcmcout.perm@S), ncol(mcmcout@S))
      checkTrue(!any(mcmcout@S == mcmcout.perm@S), "check40")
      checkEquals(nrow(mcmcout.perm@NK), nrow(mcmcout@NK))
      checkEquals(ncol(mcmcout.perm@NK), ncol(mcmcout@NK))
      checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check43")
      ## Test exception ##
      perm.index    <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
      checkException(swapElements(mcmcout, perm.index), silent = TRUE)
}

