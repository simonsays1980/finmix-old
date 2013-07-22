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
}
