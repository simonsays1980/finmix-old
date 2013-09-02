## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with finmix. If not, see <http://www.gnu.org/licenses/>.

"mcmcestimate" <- function(mcmcout, method = "kmeans", fdata = NULL, 
                           permOut = FALSE) {
    ## Check input ##
    .check.args.Mcmcestimate(mcmcout, method, fdata, permOut)
    ## Constants
    K           <- mcmcout@model@K
    M           <- mcmcout@M
    dist        <- mcmcout@model@dist
    indicmod    <- mcmcout@model@indicmod
    ranperm     <- mcmcout@ranperm

    ## If it inherits from 'mcmcoutputbase' indicators 
    ## must be simulated.
    indicfix    <- mcmcout@model@indicfix

    ## If it inherits from 'mcmcoutputperm' it has already 
    ## identified samples
    perm        <- inherits(mcmcout, what = "mcmcoutputperm")
    
    ## Posterior Mode (MAP)
    map.index   <- .map.Mcmcestimate(mcmcout)
    map         <- .extract.Mcmcestimate(mcmcout, map.index)    
    
    ## Bayesian Maximum Likelihood (BML)
    bml.index   <- .bml.Mcmcestimate(mcmcout)
    bml         <- .extract.Mcmcestimate(mcmcout, bml.index)

    ## Ergodic average (EAVG)
    eavg        <- .eavg.Mcmcestimate(mcmcout)
    
    if (indicfix) {
        ## Ergodic average is identified
        ## 'avg.id'
        ## Posterior Std. Error. 
        sdpost      <- .sdpost.Mcmcestimate(mcmcout, perm)

        .mcmcestfix(dist = dist, K = K, M = mcmcout@M, burnin = mcmcout@burnin, 
                    ranperm = mcmcout@ranperm, relabel = "none",
                    indicmod = indicmod, map = map, bml = bml, ieavg = eavg,
                    sdpost = sdpost)
    } else {
        if (ranperm) {
            ## Ergodic average is invariant
            ## 'inv'
            ## Check if already identification has been made
            if (perm) {
                if (mcmcout@Mperm > 0) {    
                    ## Use ergodic average function on 'mcmcoutputperm'
                    ## object
                    ieavg <- .eavg.Mcmcestimate(mcmcout)
                    ## Posterior Std. Error. 
                    sdpost      <- .sdpost.Mcmcestimate(mcmcout, perm)                                       
                    .mcmcestfix(dist = dist, K = K, 
                                indicmod = indicmod, M = mcmcout@Mperm, 
                                burnin = mcmcout@burnin, ranperm = mcmcout@ranperm,
                                relabel = mcmcout@relabel, map = map, bml = bml, 
                                ieavg = ieavg, sdpost = sdpost)
                } else {
                    warning(paste("No identification possible. Not a single ",
                                  "draw is a permutation", sep = ""))
                }
            } else {
                ## Use function 'mcmcpermute' to permute the sample
                mcmcoutperm <- mcmcpermute(mcmcout, method = method, fdata = fdata)
                perm    <- TRUE
                if (mcmcoutperm@Mperm > 0) {
                    ## Use ergodic average function on 'mcmcoutputperm'
                    ## object
                    ## Build 'avg.id'
                    ieavg <- .eavg.Mcmcestimate(mcmcoutperm)
                    ## Posterior Std. Error
                    sdpost  <- .sdpost.Mcmcestimate(mcmcoutperm, perm)
                    mcmcest <- .mcmcestind(dist = dist, K = K, 
                                           indicmod = indicmod, M = mcmcoutperm@Mperm, 
                                           burnin = mcmcout@burnin, ranperm = mcmcout@ranperm, 
                                           relabel = method, map = map, bml = bml, 
                                           ieavg = ieavg, eavg = eavg, sdpost = sdpost)
                    if (permOut) {
                        return.list <- list(mcmcest = mcmcest, 
                                            mcmcoutputperm = mcmcoutperm)
                        return(return.list)
                    } else {
                        return(mcmcest)
                    }
                } else {
                    warning(paste("No identification possible. Not a single ",
                                  "draw is a permutation", sep = ""))
                }
            }
        } else { 
            ## 'eavg'
            ## Posterior Std. Error
            sdpost  <- .sdpost.Mcmcestimate(mcmcout, perm)
            .mcmcestfix(dist = dist, K = K, indicmod = indicmod,
                        M = mcmcout@M, burnin = mcmcout@burnin, 
                        ranperm = mcmcout@ranperm, relabel = "none",
                        map = map, bml = bml, ieavg = eavg, sdpost = sdpost)
        }
    }
    ## New 'mcmcestimate' object.
    
    ## In case the permOut = TRUE the mcmcout object is 
    ## returned as well in a list
}

### Private functions
### These functions are not exported.

### Checking
### Check arguments: The 'mcmcout' object must inherit from
### 'mcmcoutput' or 'mcmcoutputperm'. Argument 2 must match one
### of three permutation algorithms in 'mcmcpermute()'.
### Argument 3 must be of type logical. If any case is not true 
### an error is thrown.
".check.args.Mcmcestimate" <- function(obj, arg2, arg3, arg4)
{
    if (!inherits(obj, c("mcmcoutput", "mcmcoutputperm"))) {
        stop(paste("Wrong argument: Argument 1 must be an object ",
                   "either of class 'mcmcoutput' or of type ",
                   "'mcmcoutputperm'.", sep = ""))
    } 
    match.arg(arg2, c("kmeans", "Stephens1997a", "Stephens1997b"))
    if (!inherits(arg3, "fdata") && !is.null(arg3)) {
        stop(paste("Wrong argument: Argument 3 must be an object ",
                   "of class 'fdata'.", sep = ""))
    }
    if (!is.logical(arg4)) {
        stop("Wrong argument: Argument 4 must be of type 'logical'.")
    }
}

".map.Mcmcestimate" <- function(obj) {
    ## Take the value with the highest posterior log
    ## likelihood
    mixpost <- obj@log$mixlik + obj@log$mixprior
    mixpost.sort <- sort.int(mixpost, index.return = TRUE)
    map.index <- tail(mixpost.sort$ix, 1)
    return(as.integer(map.index))
}

".bml.Mcmcestimate" <- function(obj) {
    ## Take the value with the highest log likelihood 
    mixlik <- obj@log$mixlik
    mixlik.sort <- sort.int(mixlik, index.return = TRUE)
    bml.index <- tail(mixlik.sort$ix, 1)
    return(bml.index)
}
    
".extract.Mcmcestimate" <- function(obj, m) {
    ## Extract the 'm'th row in each slot of an mcmcout
    ## object
    K           <- obj@model@K
    dist        <- obj@model@dist
    indicfix    <- !inherits(obj, what = "mcmcoutputbase")
    if (dist == "poisson") {
        par.est <- list(lambda = as.array(obj@par$lambda[m, ]))
    }
    if(!indicfix && K > 1) {
        weight.est  <- as.array(obj@weight[m, ])
        est.list    <- list(par = par.est,  weight = weight.est, 
                            log = obj@log$mixlik[m])
        return(est.list)
    }
    est.list <- list(par = par.est, log = obj@log$mixlik[m])
    return(est.list)
}

".eavg.Mcmcestimate" <- function(obj) {
    ## Check arguments ##
    dist <- obj@model@dist
    indicfix <- !inherits(obj, what = "mcmcoutputbase")
    perm <- inherits(obj, what = "mcmcoutputperm")
    if (dist == "poisson") {
        if(!perm) {
            par.eavg <- list(lambda = as.array(apply(obj@par$lambda, 
                                            2, mean, na.rm = TRUE)))
        } else {
            par.eavg <- list(lambda = as.array(apply(obj@parperm$lambda,
                                            2, mean, na.rm = TRUE)))
        }
    }
    if (indicfix) {
        eavg.list <- list(par = par.eavg)
        return(eavg.list)
    } else {
        if (perm) {
             weight.eavg <- as.array(apply(obj@weightperm,
                                           2, mean, na.rm = TRUE))
             eavg.list <- list(par = par.eavg, weight = weight.eavg)
             return(eavg.list)
        } else {
            weight.eavg = as.array(apply(obj@weight, 2, mean, 
                                         na.rm = TRUE))
            eavg.list <- list(par = par.eavg, weight = weight.eavg)
            return(eavg.list)
        }
    }
}

".sdpost.Mcmcestimate" <- function(obj, perm) 
{
    .sdpost.poisson.Mcmcestimate(obj, perm)
}

".sdpost.poisson.Mcmcestimate" <- function(obj, perm)
{
    if (perm) {
        sdpar       <- apply(obj@parperm$lambda, 2, sd)
        sdparpre    <- apply(obj@par$lambda, 2, sd)
        sdweight    <- apply(obj@weightperm, 2, sd)
        sdweightpre <- apply(obj@weight, 2, sd)
        identified      <- list(par = list(lambda = sdpar), 
                                weight = sdweight)
        unidentified    <- list(par = list(lambda = sdparpre), 
                                weight = sdweightpre)
        sdlist          <- list(identified = identified, 
                                unidentified = unidentified)    
    } else {
        sdpar       <- apply(obj@par$lambda, 2, sd)
        sdweight    <- apply(obj@weight)
        identfied   <- list(par = list(lambda = sdpar), weight = sdweight)
        sdlist      <- list(identified = identified)
    }
    return(sdlist)
}
