### --- Test Setup --- ###

if (TRUE) {
    ## Not really needed, but can be handy
    ## when writing tests
    library("RUnit")
    library("finmix")
}


"test.mcmcpm" <- function() {
    ## Set up the test
    set.seed(0)
    data <- .setUp.data(withInd = FALSE)
    model <- .setUp.model()#
    setK(model) <- 2
    setIndicFix(model) <- TRUE
    prior <- prior(hier = FALSE)
    prior <- priordefine(data, model, varargin = prior)
    mcmc <- .setUp.mcmc()
    mcmcout <- mixturemcmc(data, model, prior, mcmc)
    pm.index <- mcmc.pm(mcmcout) 
    ## Test cases ##
    checkTrue(!is.null(pm.index), "check1")
    checkTrue(is.integer(pm.index), "check2")
    checkTrue(length(pm.index) == 1, "check3")
}
