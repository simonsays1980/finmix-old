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

"mcmcpermute" <- function(mcmcout) {
    ## Check arguments ##
    .check.args.Mcmcpermute(mcmcout)
    if (mcmcout@model@K == 1) {
        return(mcmcout)
    }
    mcmcout <- .mcmcpermute.Coerce(mcmcout)
    ## Constants ##
    K           <- mcmcout@model@K
    M           <- mcmcout@M
    dist        <- mcmcout@model@dist
    indicmod    <- mcmcout@model@indicmod
   
    ## Calculate maximum a posterior estimate (MAP)
    map.index   <- mcmc.map(mcmcout)
    map         <- mcmc.extract(mcmcout, map.index)
    
    if (dist == "poisson") {
        clust.par       <- sqrt(mcmcout@par$lambda)
        clust.par       <- as.vector(clust.par)
        clust.center    <- sqrt(map$par$lambda)
    }

    ## Apply unsupervised k-means clustering to parameters
    result.clust    <- kmeans(clust.par, centers = as.vector(clust.center))
    perm.index      <- array(result.clust$clust, dim = c(M, K))
    comp.index      <- as.array(matrix(seq(1:K), nrow = M, ncol = K, 
                                       byrow = TRUE))
    keep.index      <- (t(apply(perm.index, 1, sort, FALSE)) 
                        == comp.index)
    is.perm         <- matrix(apply(comp.index, 1, all))
    nonperm         <- sum(!is.perm)
    if (nonperm < M) {
        ## Create a subsequence of the MCMC output
        mcmcout.subseq <- subseq(mcmcout, is.perm)

        ## Apply permutation suggested by kmeans clustering
        mcmcout.swap <- swapElements(mcmcout.subseq, perm.index[is.perm,])

        ## Create 'mcmcoutputperm' objects ##
        .process.Output(mcmcout, mcmcout.swap)
    } else {
        .process.Output.Empty(mcmcout)
    }
}

### Private functions.
### These functions are not exported.

### Checking
### Check arguments: Arguments must inherit either from the 
### 'mcmcoutput' class or from the 'mcmcoutputperm' class.
"check.args.Mcmcpermute" <- function(obj)
{
    if (!inherits(obj, c("mcmcoutput", "mcmcoutputperm"))) {
        stop(paste("Argument 'mcmcout' must inherit either from",
                   "class 'mcmcoutput' or from class 'mcmcoutputperm'.",
                   sep = ""))
    }
}

".mcmcpermute.Coerce" <- function(obj)
{
    ## If object is of class 'mcmcoutputperm' coerce it 
    ## to an object of class 'mcmcoutput' 
    if(inherits(obj, "mcmcoutputperm")) {
        if (class(obj) == "mcmcoutputpermfix") {
            obj <- as(obj, "mcmcoutputfix")
        } else if (class(obj) == "mcmcoutputpermfixhier") {
            obj <- as(obj, "mcmcoutputfixhier") 
        } else if (class(obj) == "mcmcoutputpermfixpost") {
            obj <- as(obj, "mcmcoutputfixpost")
        } else if (class(obj) == "mcmcoutputpermfixhierpost") {
            obj <- as(obj, "mcmcoutputfixhierpost")
        } else if (class(obj) == "mcmcoutputpermbase") {
            obj <- as(obj, "mcmcoutputbase")
        } else if (class(obj) == "mcmcoutputpermhier") {
            obj <- as(obj, "mcmcoutputhier")
        } else if (class(obj) == "mcmcoutputpermpost") {
            obj <- as(obj, "mcmcoutputpost")
        } else {
            obj <- as(obj, "mcmcoutputhierpost")
        }
    }
    return(obj)
}

".process.Output" <- function(obj, obj.swap)
{
    ## Create 'mcmcoutputperm' objects ##
    if (class(obj) == "mcmcoutputfix") {
        .mcmcoutputpermfix(obj, 
                           Mperm        = obj.swap@M,
                           parperm      = obj.swap@par,
                           logperm      = obj.swap@log)
    } else if (class(obj) == "mcmcoutputfixhier") {
        .mcmcoutputpermfixhier(obj, 
                               Mperm        = obj.swap@M,
                               parperm      = obj.swap@par,
                               logperm      = obj.swap@log)
    } else if (class(obj) == "mcmcoutputfixpost") {
        .mcmcoutputpermfixpost(obj,
                               Mperm        = obj.swap@M,
                               parperm      = obj.swap@par,
                               logperm      = obj.swap@log,
                               postperm     = obj.swap@post)
    } else if (class(obj) == "mcmcoutputfixhierpost") {
        .mcmcoutputpermfixhierpost(obj,
                                   Mperm        = obj.swap@M,
                                   parperm      = obj.swap@par,
                                   logperm      = obj.swap@log,
                                   postperm     = obj.swap@post)
    } else if (class(obj) == "mcmcoutputbase") {
        .mcmcoutputpermbase(obj,
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,                            
                            weightperm   = obj.swap@weight,
                            logperm      = obj.swap@log,
                            entropyperm  = obj.swap@entropy,
                            STperm       = obj.swap@ST,
                            Sperm        = obj.swap@S,
                            NKperm       = obj.swap@NK)          
    } else if (class(obj) == "mcmcoutputhier") {
        .mcmcoutputpermhier(obj,
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,
                            weightperm   = obj.swap@weight,
                            logperm      = obj.swap@log,
                            entropyperm  = obj.swap@entropy,
                            STperm       = obj.swap@ST,
                            Sperm        = obj.swap@S,
                            NKperm       = obj.swap@NK)
    } else if (class(obj) == "mcmcoutputpost") {
        .mcmcoutputpermpost(obj,
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,
                            weightperm   = obj.swap@weight,
                            logperm      = obj.swap@log,
                            postperm     = obj.swap@post,
                            entropyperm  = obj.swap@entropy,
                            STperm       = obj.swap@ST,
                            Sperm        = obj.swap@S,
                            NKperm       = obj.swap@NK)           
    } else {
        .mcmcoutputpermhierpost(obj,
                                Mperm        = obj.swap@M,
                                parperm      = obj.swap@par,
                                weightperm   = obj.swap@weight,
                                logperm      = obj.swap@log,
                                postperm     = obj.swap@post,
                                entropyperm  = obj.swap@entropy,
                                STperm       = obj.swap@ST,
                                Sperm        = obj.swap@S,
                                NKperm       = obj.swap@NK)           
    }
}

".process.Output.Empty" <- function(obj)
{
    warning(paste("Not a single draw is a permutation in the",
                  "function 'mcmcpermute()'.", sep = ""))
    ## Create 'mcmcoutputperm' objects ##
    if (class(obj) == "mcmcoutputfix") {
        .mcmcoutputpermfix(obj) 
    } else if (class(obj) == "mcmcoutputfixhier") {
        .mcmcoutputpermfixhier(obj) 
    } else if (class(obj) == "mcmcoutputfixpost") {
        .mcmcoutputpermfixpost(obj)
    } else if (class(obj) == "mcmcoutputfixhierpost") {
        .mcmcoutputpermfixhierpost(obj)
    } else if (class(obj) == "mcmcoutputbase") {
        .mcmcoutputpermbase(obj)
    } else if (class(obj) == "mcmcoutputhier") {
        .mcmcoutputpermhier(obj,)
    } else if (class(obj) == "mcmcoutputpost") {
        .mcmcoutputpermpost(obj) 
    } else {
        .mcmcoutputpermhierpost(obj)
    }
}
