"mcmcestimate" <- function(mcmcout, returnOut = FALSE) {
    ## Check input ##
    if (!inherits(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must be either of type 
             'mcmcoutput' or of type 'mcmcoutputperm'.")
    } 
    if (!is.logical(returnOut)) {
        stop("Argument 'returnOut' must be of type 'logical'.")
    }

    ## Constants
    K           <- mcmcout@model@K
    M           <- mcmcout@M
    dist        <- mcmcout@model@dist
    indicmod    <- mcmcout@model@indicmod
    ranperm     <- mcmcout@ranperm

    ## If it inherits from 'mcmcoutputbase' indicators 
    ## must be simulated.
    indicfix    <- !inherits(mcmcout, what = "mcmcoutputbase")

    ## If it inherits from 'mcmcoutputperm' it has already 
    ## identified samples
    perm        <- inherits(mcmcout, what = "mcmcoutputperm")
    
    ## Posterior Mode (MAP)
    map.index   <- mcmc.map(mcmcout)
    map         <- mcmc.extract(mcmcout, map.index) 
    
    ## Bayesian Maximum Likelihood (BML)
    bml.index   <- mcmc.bml(mcmcout)
    bml         <- mcmc.extract(mcmcout, bml.index)

    ## Ergodic average (EAVG)
    eavg        <- mcmc.eavg(mcmcout)

    if (indicfix) {
        ## Ergodic average is identified
        ## 'avg.id'
        mcmcest <- new("mcmcestfix", dist = dist, K = K, indicmod = indicmod,
                       map = map, bml = bml, ieavg = eavg)
        return(mcmcest)
    } else {
        if (ranperm) {
            ## Ergodic average is invariant
            ## 'inv'
            ## Check if already identification has been made
            if (perm) {
                if (mcmcout@Mperm > 0) {    
                    ## Use ergodic average function on 'mcmcoutputperm'
                    ## object
                    ieavg <- mcmc.eavg(mcmcout) 
                    mcmcest <- new("mcmcestfix", dist = dist, K = K, 
                                   indicmod = indicmod, map = map, bml = bml,
                                   ieavg = ieavg)
                    return(mcmcest)
                } else {
                    warning("No identification possible. Not a single draw in
                         'mcmcoutputperm' object is a permutation")
                }
            } else {
                ## Use function 'mcmcpermute' to permute the sample
                mcmcoutperm <- mcmcpermute(mcmcout)
                if (mcmcoutperm@Mperm > 0) {
                    ## Use ergodic average function on 'mcmcoutputperm'
                    ## object
                    ## Build 'avg.id'
                    ieavg <- mcmc.eavg(mcmcout)
                    mcmcest <- new("mcmcestind", dist = dist, K = K, 
                                   indicmod = indicmod, map = map, bml = bml,
                                   ieavg = ieavg, eavg = eavg)
                    if (returnOut) {
                        return.list <- list(mcmcest = mcmcest, 
                                            mcmcoutputperm = mcmcoutperm)
                    } else {
                        return(mcmcest)
                    }
                } else {
                    warning("No identification possible. Not a single draw in
                         'mcmcoutputperm' object is a permutation.")
                }
            }
        } else { 
            ## 'eavg'
            mcmcest <- new("mcmcestfix", dist = dist, K = K, indicmode = indicmod,
                          pm = pm, bml = bml, ieavg = eavg)
            return(mcmcest)
        }
    }
    ## New 'mcmcestimate' object.
    
    ## In case the returnOut = TRUE the mcmcout object is 
    ## returned as well in a list
}

"mcmc.map" <- function(mcmcout) {
    ## Take the value with the highest posterior log
    ## likelihood
    ## Check arguments ##
    if (!inherits(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must be either of type
             'mcmcoutput' or of type 'mcmcoutputperm'.")
    }
    mixpost <- mcmcout@log$mixlik + mcmcout@log$mixprior
    mixpost.sort <- sort.int(mixpost, index.return = TRUE)
    map.index <- tail(mixpost.sort$ix, 1)
    return(as.integer(map.index))
}

"mcmc.bml" <- function(mcmcout) {
    ## Take the value with the highest log likelihood 
    ## Check argument ##
    if (!inherits(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must be either of type 
             'mcmcoutput' or of type 'mcmcoutputperm'.")
    }
    mixlik <- mcmcout@log$mixlik
    mixlik.sort <- sort.int(mixlik, index.return = TRUE)
    bml.index <- tail(mixlik.sort$ix, 1)
    return(bml.index)
}
    
"mcmc.extract" <- function(mcmcout, m) {
    ## Extract the 'm'th row in each slot of an mcmcout
    ## object
    ## Check arguments ##
    if (!inherits(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must be either of type
             'mcmcoutput' or of type 'mcmcoutputperm'.")
    }
    K <- mcmcout@model@K
    dist <- mcmcout@model@dist
    indicfix <- !inherits(mcmcout, what = "mcmcoutputbase")
    if (dist == "poisson") {
        par.est <- list(lambda = as.array(mcmcout@par$lambda[m, ], 
                                          dim = c(1, K)))
    }
    if(!indicfix && K > 1) {
        weight.est <- mcmcout@weight[m, ]
        est.list <- list(par = par.est,  weight = weight.est)
        return(map.list)
    }
    est.list <- list(par = par.est)
    return(est.list)
}

"mcmc.eavg" <- function(mcmcout) {
    ## Check arguments ##
    if (!inherits(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must be either of type 
             'mcmcoutput' or of type 'mcmcoutputperm'.")
    }
    dist <- mcmcout@model@dist
    indicfix <- !inherits(mcmcout, what = "mcmcoutputbase")
    perm <- inherits(mcmcout, what = "mcmcoutputperm")
    if (dist == "poisson") {
        if(!perm) {
            par.eavg <- list(lambda = apply(mcmcout@par$lambda, 
                                            2, mean, na.rm = TRUE))
        } else {
            par.eavg <- list(lambda = apply(mcmcout@parperm$lambda,
                                            2, mean, na.rm = TRUE))
        }
    }
    if (indicfix) {
        eavg.list <- list(par = par.eavg)
        return(eavg.list)
    } else {
        if (perm) {
             weight.eavg <- as.array(apply(mcmcout@weightperm,
                                           2, mean, na.rm = TRUE))
             eavg.list <- list(par = par.eavg, weight = weight.eavg)
             return(eavg.list)
        } else {
            weight.eavg = as.array(apply(mcmcout@weight, 2, mean, 
                                         na.rm = TRUE))
            eavg.list <- list(par = par.eavg, weight = weight.eavg)
            return(eavg.list)
        }
    }
}
