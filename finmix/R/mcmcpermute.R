"mcmcpermute" <- function(mcmcout) {
    ## Check arguments ##
    if (!inherits(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must inherit either from class
             'mcmcoutput' or from class 'mcmcoutputperm'.")
    }
    if (mcmcout@model@K == 1) {
        return(mcmcout)
    }
    ## If object is of class 'mcmcoutputperm' coerce it 
    ## to an object of class 'mcmcoutput' 
    if(inherits(mcmcout, "mcmcoutputperm")) {
        if (class(mcmcout) == "mcmcoutputpermfix") {
            mcmcout <- as(mcmcout, "mcmcoutputfix")
        } else if (class(mcmcout) == "mcmcoutputpermfixhier") {
            mcmcout <- as(mcmcout, "mcmcoutputfixhier") 
        } else if (class(mcmcout) == "mcmcoutputpermfixpost") {
            mcmcout <- as(mcmcout, "mcmcoutputfixpost")
        } else if (class(mcmcout) == "mcmcoutputpermfixhierpost") {
            mcmcout <- as(mcmcout, "mcmcoutputfixhierpost")
        } else if (class(mcmcout) == "mcmcoutputpermbase") {
            mcmcout <- as(mcmcout, "mcmcoutputbase")
        } else if (class(mcmcout) == "mcmcoutputpermhier") {
            mcmcout <- as(mcmcout, "mcmcoutputhier")
        } else if (class(mcmcout) == "mcmcoutputpermpost") {
            mcmcout <- as(mcmcout, "mcmcoutputpost")
        } else {
            mcmcout <- as(mcmcout, "mcmcoutputhierpost")
        }
    }

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
        if (class(mcmcout) == "mcmcoutputfix") {
            mcmcoutperm <- new("mcmcoutputpermfix", mcmcout, 
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               logperm      = mcmcout.swap@log)
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputfixhier") {
            mcmcoutperm <- new("mcmcoutputpermfixhier", mcmcout, 
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               logperm      = mcmcout.swap@log)
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputfixpost") {
            mcmcoutperm <- new("mcmcoutputpermfixpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               logperm      = mcmcout.swap@log,
                               postperm     = mcmcout.swap@post)
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputfixhierpost") {
            mcmcoutperm <- new("mcmcoutputpermfixhierpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               logperm      = mcmcout.swap@log,
                               postperm     = mcmcout.swap@post)
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputbase") {
            mcmcoutperm <- new("mcmcoutputpermbase", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               logperm      = mcmcout.swap@log,
                               entropyperm  = mcmcout.swap@entropy,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputhier") {
            mcmcoutperm <- new("mcmcoutputpermhier", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               logperm      = mcmcout.swap@log,
                               entropyperm  = mcmcout.swap@entropy,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputpost") {
            mcmcoutperm <- new("mcmcoutputpermpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               logperm      = mcmcout.swap@log,
                               postperm     = mcmcout.swap@post,
                               entropyperm  = mcmcout.swap@entropy,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        } else {
            mcmcoutperm <- new("mcmcoutputpermhierpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               logperm      = mcmcout.swap@log,
                               postperm     = mcmcout.swap@post,
                               entropyperm  = mcmcout.swap@entropy,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        }
    } else {
        Mperm <- 0
        warning("Not a single draw is a permutation in the function 
                'mcmcpermute()'.")
         ## Create 'mcmcoutputperm' objects ##
        if (class(mcmcout) == "mcmcoutputfix") {
            mcmcoutperm <- new("mcmcoutputpermfix", mcmcout, 
                               Mperm        = Mperm,
                               parperm      = list(),
                               logperm      = list())
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputfixhier") {
            mcmcoutperm <- new("mcmcoutputpermfixhier", mcmcout, 
                               Mperm        = Mperm,
                               parperm      = list(),
                               logperm      = list())
            return(mcmcoutputperm)
        } else if (class(mcmcout) == "mcmcoutputfixpost") {
            mcmcoutperm <- new("mcmcoutputpermfixpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               logperm      = list(),
                               postperm     = list())
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputfixhierpost") {
            mcmcoutperm <- new("mcmcoutputpermfixhierpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               logperm      = list(),
                               postperm     = list())
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputbase") {
            mcmcoutperm <- new("mcmcoutputpermbase", mcmcout,
                               Mperm        = MPerm,
                               parperm      = list(),
                               weightperm   = array(),
                               logperm      = list(),
                               entropyperm  = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputhier") {
            mcmcoutperm <- new("mcmcoutputpermhier", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               weightperm   = array(),
                               logperm      = list(),
                               entropyperm  = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        } else if (class(mcmcout) == "mcmcoutputpost") {
            mcmcoutperm <- new("mcmcoutputpermpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               weightperm   = array(),
                               logperm      = list(),
                               postperm     = list(),
                               entropyperm  = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        } else {
            mcmcoutperm <- new("mcmcoutputpermhierpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               weightperm   = array(),
                               logperm      = list(),
                               postperm     = list(),
                               entropyperm  = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        }
       
    }
}



