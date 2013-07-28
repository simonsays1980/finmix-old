### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

".setUp.data" <- function(withInd = FALSE) {
        ## Get path ##
        pkg <- "finmix"
        if (Sys.getenv("RCMDCHECK") == FALSE) {
            data.path <- file.path(getwd(), "..", 
                                   "data", "poisson.data.csv")        
        } else {
            data.path <- system.file(package = pkg, 
                                     'data/poisson.data.csv')
        }
        data <- read.csv(data.path, header = FALSE, sep = ",")
        if(withInd) {
            if (Sys.getenv("RCMDCHECK") == FALSE) {

                ind.path <- file.path(getwd(), "..", 
                                       "data", 
                                       "poisson.ind.csv")        
            } else {
                ind.path <- system.file(package = pkg, 
                                        'data/poisson.ind.csv')
            }               
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check11")
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check18")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout@hyper$b[1:5] == mcmcout.sub@hyper$b), "check21")
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputfixpost.poisson" <- function() {
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
      sub.index             <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, 
                                      ncol = 1) 
      mcmcout.sub           <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout.sub@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@post$par$a), 5)
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5]), "check7")
      checkEquals(nrow(mcmcout.sub@post$par$b), 5)
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5]), "check10")
       
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check14")
      checkEquals(nrow(mcmcout.sub@post$par$a), nrow(mcmcout@post$par$a[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check17")
      checkEquals(nrow(mcmcout.sub@post$par$b), nrow(mcmcout@post$par$b[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check20")
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check24")
      checkEquals(nrow(mcmcout.sub@post$par$a), nrow(mcmcout@post$par$a[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check27")
      checkEquals(nrow(mcmcout.sub@post$par$b), nrow(mcmcout@post$par$b[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check30")
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputfixhierpost.poisson" <- function() {
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
      sub.index             <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, 
                                      ncol = 1) 
      mcmcout.sub           <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout.sub@par$lambda), ncol(mcmcout@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@post$par$a), 5)
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5]))
      checkEquals(nrow(mcmcout.sub@post$par$b), 5)
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5]), "check10")
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check17")
      checkEquals(nrow(mcmcout.sub@post$par$a), nrow(mcmcout@post$par$a[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check20")
      checkEquals(nrow(mcmcout.sub@post$par$b), nrow(mcmcout@post$par$b[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check23")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout@hyper$b[1:5] == mcmcout.sub@hyper$b), "check26")
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check30")
      checkEquals(nrow(mcmcout.sub@post$par$a), nrow(mcmcout@post$par$a[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check33")
      checkEquals(nrow(mcmcout.sub@post$par$b), nrow(mcmcout@post$par$b[1:5,]))
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check36")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout@hyper$b[1:5] == mcmcout.sub@hyper$b), "check39")
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputbase.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data          <- .setUp.data(withInd = TRUE)
      model         <- .setUp.model()
      setK(model)   <- 2
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmc          <- .setUp.mcmc()
      prior         <- prior(hier = FALSE)
      prior         <- priordefine(data, model, varargin = prior)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check7")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check10")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check13")
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check17")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check20")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check23")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check26")
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputhier.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data          <- .setUp.data(withInd = TRUE)
      model         <- .setUp.model()
      setK(model)   <- 2
      prior         <- priordefine(data, model)
      mcmc          <- .setUp.mcmc()
      prior         <- priordefine(data, model)
      mcmcout       <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check7")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check10")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check13")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout.sub@hyper$b == mcmcout@hyper$b[1:5]))
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check17")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check20")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check23")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check26")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout.sub@hyper$b == mcmcout@hyper$b[1:5]))
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputpost.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 2
      prior                 <- prior(hier = FALSE)
      prior                 <- priordefine(data, model, varargin = prior)
      mcmc                  <- .setUp.mcmc()
      setStorepost(mcmc)    <- TRUE
      prior                 <- priordefine(data, model)
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check7")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check10")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check13")
      checkEquals(nrow(mcmcout.sub@post$par$a), 5)
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check16")
      checkEquals(nrow(mcmcout.sub@post$par$b), 5)
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check19")
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check243")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check26")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check29")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check32")
      checkEquals(nrow(mcmcout.sub@post$par$a), 5)
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check35")
      checkEquals(nrow(mcmcout.sub@post$par$b), 5)
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check38")
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

"test.subseq.mcmcoutputhierpost.poisson" <- function() {
      ## --- Check for K = 2 --- ##
      set.seed(0)
      data                  <- .setUp.data(withInd = TRUE)
      model                 <- .setUp.model()
      setK(model)           <- 2
      prior                 <- priordefine(data, model)
      mcmc                  <- .setUp.mcmc()
      setStorepost(mcmc)    <- TRUE
      prior                 <- priordefine(data, model)
      mcmcout               <- mixturemcmc(data, model, prior, mcmc)
      sub.index     <- matrix(seq(1:mcmc@M) < 6, nrow = mcmc@M, ncol = 1)
      mcmcout.sub   <- subseq(mcmcout, sub.index)
      ## Test cases ##
      checkEquals(mcmcout.sub@M, 5)
      checkEquals(nrow(mcmcout.sub@par$lambda), 5)
      checkEquals(ncol(mcmcout@par$lambda), ncol(mcmcout.sub@par$lambda))
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check4")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check7")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check10")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check13")
      checkEquals(nrow(mcmcout.sub@post$par$a), 5)
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check16")
      checkEquals(nrow(mcmcout.sub@post$par$b), 5)
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check19")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout.sub@hyper$b == mcmcout@hyper$b[1:5]))
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
      checkTrue(all(mcmcout@par$lambda[1:5,] == mcmcout.sub@par$lambda), "check23")
      checkEquals(nrow(mcmcout.sub@weight), 5)
      checkEquals(ncol(mcmcout.sub@weight), ncol(mcmcout@weight))
      checkTrue(all(mcmcout@weight[1:5,] == mcmcout.sub@weight), "check26")
      checkEquals(nrow(mcmcout.sub@ST), 5)
      checkEquals(ncol(mcmcout.sub@ST), 1)
      checkTrue(all(mcmcout.sub@ST == mcmcout@ST[1:5]), "check29")
      checkEquals(nrow(mcmcout.sub@NK), 5)
      checkEquals(ncol(mcmcout.sub@NK), ncol(mcmcout@NK))
      checkTrue(all(mcmcout.sub@NK == mcmcout@NK[1:5,]), "check32")
      checkEquals(nrow(mcmcout.sub@post$par$a), 5)
      checkEquals(ncol(mcmcout.sub@post$par$a), ncol(mcmcout@post$par$a))
      checkTrue(all(mcmcout.sub@post$par$a == mcmcout@post$par$a[1:5,]), "check35")
      checkEquals(nrow(mcmcout.sub@post$par$b), 5)
      checkEquals(ncol(mcmcout.sub@post$par$b), ncol(mcmcout@post$par$b))
      checkTrue(all(mcmcout.sub@post$par$b == mcmcout@post$par$b[1:5,]), "check38")
      checkEquals(nrow(mcmcout.sub@hyper$b), 5)
      checkEquals(ncol(mcmcout.sub@hyper$b), ncol(mcmcout@hyper$b))
      checkTrue(all(mcmcout.sub@hyper$b == mcmcout@hyper$b[1:5]))
      ## Test exception ##
      sub.index    <- matrix(seq(1:50) < 6, nrow = 50, ncol = model@K)
      checkException(subseq(mcmcout, sub.index), silent = TRUE)
}

