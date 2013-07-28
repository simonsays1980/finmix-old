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

"test.mcmcpermute.mcmcoutputfix.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()
    setK(model) <- 1
    setIndicfix(model) <- TRUE
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfix"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model) <- 2
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfix"), "check2")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check3")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check4")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check9")
    ## --- Check K = 3 --- ##
    setK(model) <- 3
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfix"), "check10")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check11")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check12")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check13")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check17")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfix"), "check18")
}

"test.mcmcpermute.mcmcoutputfixhier.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()
    setK(model) <- 1
    setIndicfix(model) <- TRUE
    prior <- priordefine(data, model)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfixhier"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model) <- 2
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhier"), "check2")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check3")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check4")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check9")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhier"), "check10")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check11")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check12")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check13")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check17")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhier"), "check18")
}

"test.mcmcpermute.mcmcoutputfixpost.poisson" <- function() {
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
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfixpost"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model)     <- 2
    prior           <- prior(hier = FALSE)
    prior           <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixpost"), "check2")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check3")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check4")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check5")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check6")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check10")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model)     <- 3
    prior           <- prior(hier = FALSE)
    prior           <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixpost"), "check13")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check14")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check15")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check16")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check17")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check21")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixpost"), "check24")
}

"test.mcmcpermute.mcmcoutputfixhierpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data                <- .setUp.data(withInd = TRUE)
    model               <- .setUp.model()
    setK(model)         <- 1
    setIndicfix(model)  <- TRUE
    prior               <- priordefine(data, model)
    mcmc                <- .setUp.mcmc()
    setStorepost(mcmc)  <- TRUE
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfixhierpost"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model)     <- 2
    prior           <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhierpost"), "check2")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check3")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check4")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check5")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check6")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check10")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model)     <- 3
    prior           <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhierpost"), "check13")
    checkTrue("Mperm" %in% slotNames(mcmcoutperm), "check14")
    checkTrue("parperm" %in% slotNames(mcmcoutperm), "check15")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check16")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check17")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkTrue(is.integer(mcmcoutperm@Mperm), "check21")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhierpost"), "check24")
}

"test.mcmcpermute.mcmcoutputbase.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()
    setK(model) <- 1
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfix"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model) <- 2
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermbase"))
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check7")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check12")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl6")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check17")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check19")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check20")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check22")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermbase"), "check23")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check25")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check26")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check28")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check29")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check33")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check35")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check37")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check38")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check40")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check41")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check43")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermbase"), "check44")
}

"test.mcmcpermute.mcmcoutputhier.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()
    setK(model) <- 1
    prior <- priordefine(data, model)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfixhier"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model) <- 2
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhier"))
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check7")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check12")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl6")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check17")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check19")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check20")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check22")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhier"), "check23")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check25")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check26")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check28")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check29")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check33")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check35")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check37")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check38")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check40")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check41")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check43")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhier"), "check44")
}

"test.mcmcpermute.mcmcoutputpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()
    setK(model) <- 1
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    setStorepost(mcmc) <- TRUE
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfixpost"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model) <- 2
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermpost"))
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check7")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check12")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl6")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check17")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check19")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check20")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check22")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check23")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    checkEquals(ncol(mcmcoutperm@postperm$weight),
                ncol(mcmcout@post$weight))
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermpost"), "check27")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check29")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check30")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check32")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check33")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check37")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check39")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check41")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check42")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check44")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check45")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check47")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check48")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    checkEquals(ncol(mcmcoutperm@postperm$weight),
                ncol(mcmcout@post$weight))
    checkException(mcmcpermute(prior), silent = TRUE)
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermpost"), "check52")
}

"test.mcmcpermute.mcmcoutputhierpost.poisson" <- function() {
    ## --- Check K = 1 --- ##
    set.seed(0)
    data <- .setUp.data(withInd = TRUE)
    model <- .setUp.model()
    setK(model) <- 1
    prior <- priordefine(data, model)
    mcmc <- .setUp.mcmc()
    setStorepost(mcmc) <- TRUE
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputfixhierpost"), "check1")
    ## --- Check K = 2 --- ##
    set.seed(0)
    setK(model) <- 2
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check2")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check7")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check12")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl6")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check17")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check19")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check20")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check22")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check23")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    checkEquals(ncol(mcmcoutperm@postperm$weight),
                ncol(mcmcout@post$weight))
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check27")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check29")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check30")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check32")
    checkTrue("logperm" %in% slotNames(mcmcoutperm), "check33")
    checkEquals(ncol(mcmcoutperm@logperm$mixlik), 1)
    checkEquals(ncol(mcmcoutperm@logperm$mixprior), 1)
    checkEquals(ncol(mcmcoutperm@logperm$cdpost), 1)
    checkTrue("entropyperm" %in% slotNames(mcmcoutperm), "check37")
    checkEquals(ncol(mcmcoutperm@entropyperm), 1)
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check39")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check41")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check42")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check44")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check45")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check47")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check48")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    checkEquals(ncol(mcmcoutperm@postperm$weight),
                ncol(mcmcout@post$weight))
    checkException(mcmcpermute(prior), silent = TRUE)
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check52")
    ## --- Check ranperm --- ##
    set.seed(0)
    setRanperm(mcmc) <- TRUE
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check42")
}

