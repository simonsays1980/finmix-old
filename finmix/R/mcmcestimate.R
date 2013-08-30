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

"mcmcestimate" <- function(mcmcout, method = "kmeans", permOut = FALSE) {
    ## Check input ##
    .check.args.Mcmcestimate(mcmcout, method, permOut)
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
                    warning(paste("No identification possible. Not a single ",
                                  "draw is a permutation", sep = ""))
                }
            } else {
                ## Use function 'mcmcpermute' to permute the sample
                mcmcoutperm <- mcmcpermute(mcmcout, method)
                if (mcmcoutperm@Mperm > 0) {
                    ## Use ergodic average function on 'mcmcoutputperm'
                    ## object
                    ## Build 'avg.id'
                    ieavg <- mcmc.eavg(mcmcoutperm)
                    mcmcest <- new("mcmcestind", dist = dist, K = K, 
                                   indicmod = indicmod, map = map, bml = bml,
                                   ieavg = ieavg, eavg = eavg)
                    if (returnOut) {
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
            mcmcest <- new("mcmcestfix", dist = dist, K = K, indicmod = indicmod,
                          map = map, bml = bml, ieavg = eavg)
            return(mcmcest)
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
".check.args.Mcmcestimate" <- function(obj, arg2, arg3)
{
    if (!inherits(obj, c("mcmcoutput", "mcmcoutputperm"))) {
        stop(paste("Wrong argument: Argument 1 must be either of type ", 
                   "'mcmcoutput' or of type 'mcmcoutputperm'.",
                   sep = ""))
    } 
    match.arg(arg2, c("kmeans", "Stephens1997a", "stephens1997b"))
    if (!is.logical(arg3)) {
        stop("Wrong argument: Argument 3 must be of type 'logical'.")
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
        est.list    <- list(par = par.est,  weight = weight.est)
        return(est.list)
    }
    est.list <- list(par = par.est)
    return(est.list)
}

".eavg.Mcmcestimate" <- function(obj) {
    ## Check arguments ##
    dist <- obj@model@dist
    indicfix <- !inherits(obj, what = "mcmcoutputbase")
    perm <- inherits(mcmcout, what = "mcmcoutputperm")
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
