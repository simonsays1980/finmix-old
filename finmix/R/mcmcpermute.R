"mcmcpermute" <- function(mcmcout) {
    ## Check arguments ##
    if (!is(mcmcout, c("mcmcoutput", "mcmcoutputperm"))) {
        stop("Argument 'mcmcout' must be either of type
             'mcmcoutput' or of type 'mcmcoutputperm'.")
    }
    ## If object is of class 'mcmcoutputperm' coerce it 
    ## to an object of class 'mcmcoutput' 
    if(is(mcmcout, "mcmcoutputperm")) {
        if (is(mcmcout, "mcmcoutputpermfix")) {
            mcmcout <- as(mcmcout, "mcmcoutputfix")
        } else if (is(mcmcout, "mcmcoutputpermfixhier")) {
            mcmcout <- as(mcmcout, "mcmcoutputfixhier") 
        } else if(is(mcmcout, "mcmcoutputpermfixpost")) {
            mcmcout <- as(mcmcout, "mcmcoutputfixpost")
        } else if (is(mcmcout, "mcmcoutputpermfixhierpost")) {
            mcmcout <- as(mcmcout, "mcmcoutputfixhierpost")
        } else if (is(mcmcout, "mcmcoutputpermbase")) {
            mcmcout <- as(mcmcout, "mcmcoutputbase")
        } else if (is(mcmcout, "mcmcoutputpermhier")) {
            mcmcout <- as(mcmcout, "mcmcoutputhier")
        } else if (is(mcmcout, "mcmcoutputpermpost")) {
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
        clust.center    <- sqrt(map$lambda)
    }

    ## Apply unsupervised k-means clustering to parameters
    result.clust    <- kmeans(clust.par, centers = clust.center)
    perm.index      <- array(result.clust$clust, dim = c(M, K))
    comp.index      <- as.array(matrix(seq(1:K), nrow = M, ncol = K, 
                                       byrow = TRUE))
    keep.index      <- (t(apply(perm.ind, 1, sort, FALSE)) 
                        == comp.index)
    is.perm         <- apply(comp, 1, all)
    nonperm         <- sum(!is.perm)
    
    if (nonperm < M) {
        ## Create a subsequence of the MCMC output
        mcmcout.subseq <- subseq(mcmcout, is.perm)

        ## Apply permutation suggested by kmeans clustering
        mcmcout.swap <- swapElements(mcmcout.subseq, perm.index[is.perm,])

        ## Create 'mcmcoutputperm' objects ##
        if (is(mcmcout, "mcmcoutputfix")) {
            mcmcoutperm <- new("mcmcoutputpermfix", mcmcout, 
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par)
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputfixhier")) {
            mcmcoutperm <- new("mcmcoutputpermfixhier", mcmcout, 
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par)
            return(mcmcoutputperm)
        } else if (is(mcmcout, "mcmcoutputfixpost")) {
            mcmcoutperm <- new("mcmcoutputpermfixpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               postperm     = mcmcout.swap@post)
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputfixhierpost")) {
            mcmcoutperm <- new("mcmcoutputpermfixhierpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               postperm     = mcmcout.swap@post)
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputbase")) {
            mcmcoutperm <- new("mcmcoutputpermbase", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputhier")) {
            mcmcoutperm <- new("mcmcoutputpermhier", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputpost")) {
            mcmcoutperm <- new("mcmcoutputpermpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               weightperm   = mcmcout.swap@weight,
                               postperm     = mcmcout.swap@post,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        } else {
            mcmcoutperm <- new("mcmcoutputpermhierpost", mcmcout,
                               Mperm        = mcmcout.swap@M,
                               parperm      = mcmcout.swap@par,
                               postperm     = mcmcout.swap@post,
                               weightperm   = mcmcout.swap@weight,
                               STperm       = mcmcout.swap@ST,
                               Sperm        = mcmcout.swap@S,
                               NKperm       = mcmcout.swap@NK)
            return(mcmcoutperm)
        }
    } else {
        Mperm <- 0
        warning("Not a single draw is a permutation in the function 
                'mcmcpermute()'")
         ## Create 'mcmcoutputperm' objects ##
        if (is(mcmcout, "mcmcoutputfix")) {
            mcmcoutperm <- new("mcmcoutputpermfix", mcmcout, 
                               Mperm        = Mperm,
                               parperm      = list())
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputfixhier")) {
            mcmcoutperm <- new("mcmcoutputpermfixhier", mcmcout, 
                               Mperm        = Mperm,
                               parperm      = list())
            return(mcmcoutputperm)
        } else if (is(mcmcout, "mcmcoutputfixpost")) {
            mcmcoutperm <- new("mcmcoutputpermfixpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               postperm     = list())
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputfixhierpost")) {
            mcmcoutperm <- new("mcmcoutputpermfixhierpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               postperm     = list())
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputbase")) {
            mcmcoutperm <- new("mcmcoutputpermbase", mcmcout,
                               Mperm        = MPerm,
                               parperm      = list(),
                               weightperm   = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputhier")) {
            mcmcoutperm <- new("mcmcoutputpermhier", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               weightperm   = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        } else if (is(mcmcout, "mcmcoutputpost")) {
            mcmcoutperm <- new("mcmcoutputpermpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               weightperm   = array(),
                               postperm     = list(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        } else {
            mcmcoutperm <- new("mcmcoutputpermhierpost", mcmcout,
                               Mperm        = Mperm,
                               parperm      = list(),
                               postperm     = list(),
                               weightperm   = array(),
                               STperm       = array(),
                               Sperm        = array(),
                               NKperm       = array())
            return(mcmcoutperm)
        }
       
    }
}



