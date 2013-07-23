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
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfix"))
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    ## --- Check K = 3 --- ##
    setK(model) <- 3
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfix"), "check5")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check7")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfix"), "check8")
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
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhier"), "check5")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check7")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhier"), "check10")
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
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
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
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixpost"), "check7")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check9")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixpost"), "check10")
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
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
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
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhierpost"), "check7")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check9")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermfixhierpost"), "check12")
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
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl0")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check11")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check13")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check16")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermbase"), "check17")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check19")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check20")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check22")
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check23")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check25")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check26")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check28")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check29")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check31")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermbase"), "check32")
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
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl0")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check11")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check13")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check16")
    ## --- Check K = 3 --- ##
    set.seed(0)
    setK(model) <- 3
    prior <- priordefine(data, model)
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhier"), "check17")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check19")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check20")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check22")
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check23")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check25")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check26")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check28")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check29")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check31")
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhier"), "check32")
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
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl0")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check11")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check13")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check16")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check17")
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
    checkTrue(is(mcmcoutperm, "mcmcoutputpermpost"), "check21")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check23")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check24")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check26")
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check27")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check29")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check30")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check32")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check33")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check35")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check36")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    checkEquals(ncol(mcmcoutperm@postperm$weight),
                ncol(mcmcout@post$weight))
    checkException(mcmcpermute(prior), silent = TRUE)
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermpost"), "check41")
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
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"))
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check4")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check5")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check7")
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check8")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "checkl0")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check11")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check13")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check14")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check16")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check17")
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
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check21")
    checkEquals(ncol(mcmcoutperm@parperm$lambda), 
                ncol(mcmcout@par$lambda))
    checkTrue(is.integer(mcmcoutperm@Mperm), "check23")
    checkTrue("weightperm" %in% slotNames(mcmcoutperm), "check24")
    checkEquals(ncol(mcmcoutperm@weightperm), 
              ncol(mcmcout@weight))
    checkTrue(nrow(mcmcoutperm@weightperm) <=
              nrow(mcmcout@weight), "check26")
    checkTrue("STperm" %in% slotNames(mcmcoutperm), "check27")
    checkEquals(ncol(mcmcoutperm@STperm),
                ncol(mcmcout@ST))
    checkTrue(nrow(mcmcoutperm@STperm) <= nrow(mcmcout@ST), "check29")
    checkTrue("Sperm" %in% slotNames(mcmcoutperm), "check30")
    checkEquals(nrow(mcmcoutperm@Sperm),
                nrow(mcmcout@S))
    checkTrue(ncol(mcmcoutperm@Sperm) <= ncol(mcmcout@S), "check32")
    checkTrue("NKperm" %in% slotNames(mcmcoutperm), "check33")
    checkEquals(ncol(mcmcoutperm@NKperm), ncol(mcmcout@NK))
    checkTrue(nrow(mcmcoutperm@NKperm) <= nrow(mcmcout@NK), "check35")
    checkTrue("postperm" %in% slotNames(mcmcoutperm), "check36")
    checkEquals(ncol(mcmcoutperm@postperm$par$a),
                ncol(mcmcout@post$par$a))
    checkEquals(ncol(mcmcoutperm@postperm$par$b),
                ncol(mcmcout@post$par$b))
    checkEquals(ncol(mcmcoutperm@postperm$weight),
                ncol(mcmcout@post$weight))
    checkException(mcmcpermute(prior), silent = TRUE)
    ## --- Check mcmcoutputperm --- ##
    mcmcoutperm <- mcmcpermute(mcmcoutperm)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check41")
    ## --- Check ranperm --- ##
    set.seed(0)
    setRanperm(mcmc) <- TRUE
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    mcmcoutperm <- mcmcpermute(mcmcout)
    checkTrue(is(mcmcoutperm, "mcmcoutputpermhierpost"), "check42")
}

