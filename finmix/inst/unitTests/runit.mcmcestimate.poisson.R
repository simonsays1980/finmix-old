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

"test.mcmcestimate.mcmcoutputfix.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    setIndicfix(model)  <- TRUE
    prior               <- prior(hier = FALSE)
    prior               <- priordefine(data, model, varargin = prior)
    mcmc                <- .setUp.mcmc()
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0) 
    setK(model) <- 3
    prior       <- prior(hier = FALSE)
    prior       <- priordefine(data, model, varargin = prior)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputfixhier.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    setIndicfix(model)  <- TRUE
    prior               <- priordefine(data, model)
    mcmc                <- .setUp.mcmc()
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior       <- priordefine(data, model)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputfixpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    setIndicfix(model)  <- TRUE
    prior               <- prior(hier = FALSE)
    prior               <- priordefine(data, model, varargin = prior)
    mcmc                <- .setUp.mcmc()
    setStorepost(mcmc)  <- TRUE
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0) 
    setK(model) <- 3
    prior       <- prior(hier = FALSE)
    prior       <- priordefine(data, model, varargin = prior)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputfixhierpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    setIndicfix(model) <- TRUE
    prior               <- priordefine(data, model)
    mcmc                <- .setUp.mcmc()
    setStorepost(mcmc)  <- TRUE
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior       <- priordefine(data, model)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputbase.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data        <- .setUp.data(withInd = TRUE)
    model       <- .setUp.model()
    setK(model) <- 1
    setIndicfix(model) <- TRUE
    prior       <- prior(hier = FALSE)
    prior       <- priordefine(data, model, varargin = prior)
    mcmc        <- .setUp.mcmc()
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0) 
    setK(model) <- 3
    prior       <- prior(hier = FALSE)
    prior       <- priordefine(data, model, varargin = prior)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputhier.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data        <- .setUp.data(withInd = TRUE)
    model       <- .setUp.model()
    setK(model) <- 1
    prior       <- priordefine(data, model)
    mcmc        <- .setUp.mcmc()
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior       <- priordefine(data, model)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    prior               <- prior(hier = FALSE)
    prior               <- priordefine(data, model, varargin = prior)
    mcmc                <- .setUp.mcmc()
    setStorepost(mcmc)  <- TRUE
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0) 
    setK(model) <- 3
    prior       <- prior(hier = FALSE)
    prior       <- priordefine(data, model, varargin = prior)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputhierpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    prior               <- priordefine(data, model)
    mcmc                <- .setUp.mcmc()
    setStorepost(mcmc)  <- TRUE
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior       <- priordefine(data, model)
    mcmcout     <- mixturemcmc(data, model, prior, mcmc)
    mcmcest     <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")   
}

"test.mcmcestimate.mcmcoutputhierpost.ranperm.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    prior               <- priordefine(data, model)
    mcmc                <- .setUp.mcmc()
    setStorepost(mcmc)  <- TRUE
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue(!("eavg" %in% slotNames(mcmcest)), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(!is.list(mcmcest.list), "check18")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model)         <- 3
    prior               <- priordefine(data, model)
    setRanperm(mcmc)    <- TRUE
    mcmcout             <- mixturemcmc(data, model, prior, mcmc)
    mcmcest             <- mcmcestimate(mcmcout)
    checkTrue(is(mcmcest, "mcmcestfix"), "check1")
    checkTrue("dist" %in% slotNames(mcmcest), "check2")
    checkTrue("K" %in% slotNames(mcmcest), "check3")
    checkEquals(mcmcest@K, model@K)
    checkTrue("indicmod" %in% slotNames(mcmcest), "check4")
    checkTrue("map" %in% slotNames(mcmcest), "check5")
    checkTrue("bml" %in% slotNames(mcmcest), "check6")
    checkTrue("ieavg" %in% slotNames(mcmcest), "check7")
    checkTrue("eavg" %in% slotNames(mcmcest), "check8")
    checkTrue("par" %in% names(mcmcest@map), "check9")
    checkTrue("lambda" %in% names(mcmcest@map$par), "check10")
    checkEquals(dim(mcmcest@map$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@bml), "check12")
    checkTrue("lambda" %in% names(mcmcest@bml$par), "check13")
    checkEquals(dim(mcmcest@bml$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@ieavg), "check15")
    checkTrue("lambda" %in% names(mcmcest@ieavg$par), "check16")
    checkEquals(dim(mcmcest@ieavg$par$lambda), mcmcest@K)
    checkTrue("par" %in% names(mcmcest@eavg), "check17")
    checkTrue("lambda" %in% names(mcmcest@eavg$par), "check18")
    checkEquals(dim(mcmcest@eavg$par$lambda), mcmcest@K)
    ## --- returnOut = TRUE --- ##
    mcmcest.list <- mcmcestimate(mcmcout, returnOut = TRUE)
    checkTrue(is.list(mcmcest.list), "check20")   
    checkTrue("mcmcoutputperm" %in% names(mcmcest.list), "check21")
    checkTrue(is(mcmcest.list$mcmcoutputperm, "mcmcoutputperm"), "check22")
}

