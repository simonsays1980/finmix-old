### --- Test Setup --- ###

if (TRUE) {
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
            data <- data(y = as.matrix(data), S = as.matrix(ind), type = "discrete",
                         r = 1, N = nrow(data), sim = TRUE,
                         bycolumn = TRUE)
            return(data)
        } else {
                data <- data(y = as.matrix(data), type = "discrete", r = 1,
                             N = nrow(data), sim = TRUE,
                             bycolumn = TRUE)
                return(data)
        }
}

".setUp.model" <- function() {
        model <- model(dist = "poisson", K = 2)
        return(model)
}

".setUp.mcmc" <- function() {
        mcmc <- mcmc(burnin = 0, M = 100, startpar = FALSE,
                        storeS = 2, storepost = FALSE,
                        ranperm = FALSE)
        return(mcmc)
}

"test.mcmcmap" <- function() {
    ## Set up the test
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()#
    setK(model) <- 2
    setIndicfix(model) <- TRUE
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    map.index <- mcmc.map(mcmcout)
    ## Test cases ##
    checkTrue(!is.null(map.index), "check1")
    checkTrue(is.integer(map.index), "check2")
    checkTrue(length(map.index) == 1, "check3")
    checkTrue(map.index > 0, "check4")
    checkTrue(map.index < mcmcout@M, "check5")
}

"test.mcmcbml" <- function() {
    ## Set up the test
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()#
    setK(model) <- 1
    setIndicfix(model) <- TRUE
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    bml.index <- mcmc.bml(mcmcout)
    ## Test cases ##
    checkTrue(!is.null(bml.index), "check1")
    checkTrue(is.integer(bml.index), "check2")
    checkTrue(length(bml.index) == 1, "check3")
    checkTrue(bml.index > 0, "check4")
    checkTrue(bml.index < mcmcout@M, "check5")
}

"test.mcmceavg" <- function() {
    ## Set up the test
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()#
    setK(model) <- 1
    setIndicfix(model) <- TRUE
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    eavg <- mcmc.eavg(mcmcout)
    ## Test cases ##
    checkTrue(is.list(eavg), "check1")
    checkTrue(length(eavg) != 0, "check2")
    checkTrue(is.list(eavg[[1]]), "check3")
    checkTrue(is.numeric(eavg[[1]][[1]]), "check4")
    for (k in 1:model@K) {
        checkTrue(eavg[[1]][[1]][k] > 0, 
                  paste("check", 4 + k, sep = ""))
    }
}

"test.mcmcextract" <- function() {
    ## Set up the test
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()#
    setK(model) <- 1
    setIndicfix(model) <- TRUE
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    map.index <- mcmc.map(mcmcout)
    extract <- mcmc.extract(mcmcout, map.index)
    ## Test cases ##
    checkTrue(is.list(extract), "check1")
    checkTrue(length(extract) != 0, "check2")
    checkTrue(is.list(extract[[1]]), "check3")
    checkTrue(is.numeric(extract[[1]][[1]]), "check4")
    for (k in 1:model@K) {
        checkTrue(extract[[1]][k] > 0, 
                  paste("check", 4 + k, sep = ""))
    }
}
