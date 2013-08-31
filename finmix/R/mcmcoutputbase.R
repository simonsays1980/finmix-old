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

.mcmcoutputbase <- setClass("mcmcoutputbase", 
                            representation(weight 	= "array",
                                           entropy	= "array",
                                           ST 	    = "array",
                                           S 	    = "array",
                                           NK 	    = "array",
                                           clust 	= "array"),
                            contains = c("mcmcoutputfix"),
                            validity = function(object) 
                            {
                                ## else: OK
                                TRUE
                            },
                            prototype(weight    = array(),
                                      entropy   = array(),
                                      ST        = array(),
                                      S         = array(),
                                      NK        = array(),
                                      clust     = array()
                                      )
)

setMethod("show", "mcmcoutputbase", 
          function(object)
          {
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              if (!all(is.na(object@S))) {
                  cat("     S           :", 
                      paste(dim(object@S), collapse = "x"), "\n")
              }
              cat("     NK          :",
                  paste(dim(object@NK), collapse = "x"), "\n")
              cat("     clust       :",
                  paste(dim(object@clust), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plot", signature(x = "mcmcoutputbase", 
                            y = "ANY"), 
          function(x, y = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .traces.Poisson.Base(x, y)
              }
              ## log ##
              .traces.Log.Base(x, y)
          }
)

setMethod("plotHist", signature(x   = "mcmcoutputbase", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .hist.Poisson.Base(x, dev)
              }
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputbase",
                                dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if (x@model@dist == "poisson") {
                  .dens.Poisson.Base(x, dev)
              }
          }
)

setMethod("subseq", signature(object = "mcmcoutputbase", 
                              index = "array"), 
          function(object, index) 
          {
              ## Call 'subseq()' method from 'mcmcoutputfix'
              as(object, "mcmcoutputfix") <- callNextMethod(object, index)
              ## Change owned slots ##
              .subseq.Base(object, index)
          }
)

setMethod("swapElements", signature(object = "mcmcoutputbase", 
                                    index = "array"),
          function(object, index) 
          {              
              if (object@model@K == 1) {
                  return(object)
              } else {
                  ## Call method 'swapElements()' from 'mcmcoutputfix' 
                  as(object, "mcmcoutputfix") <- callNextMethod(object, index)
                  .swapElements.Base(object, index)
              }
          }
)

setMethod("getWeight", "mcmcoutputbase", 
          function(object) 
          {
              return(object@weight)
          }
)

setMethod("getEntropy", "mcmcoutputbase", 
          function(object) 
          {
              return(object@entropy)	
          }
)

setMethod("getST", "mcmcoutputbase", 
          function(object) 
          {
              return(object@ST)	
          }
)

setMethod("getS", "mcmcoutputbase", 
          function(object) 
          {
              return(object@S)	
          }
)

setMethod("getNK", "mcmcoutputbase", 
          function(object) 
          {
              return(object@NK)	
          }
)

setMethod("getClust", "mcmcoutputbase", 
          function(object) 
          {
              return(object@clust)	
          }
)

## No setters as users are not intended to manipulate ##
## this object. ##

### Private functions.
### These functions are not exported.

### Plot
### Plot traces
### Plot traces Poisson: Plots the traces for the sampled 
### Poisson parameters and the weights.
".traces.Poisson.Base" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K * 2 - 1
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    lambda <- x@par$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")                      
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]),
              cex = 0.6, line = 3)
    }
    weight <- x@weight
    for (k in 1:(K - 1)) {
        plot(weight[, k], type = "l", axes = F, 
             col = "gray47", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(eta[k = .(k)]),
              cex = 0.6, line = 3)
    }
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)                  
}

### Plot traces log-likelihood: Plots the traces of the 
### sampled log-likelihoods.
".traces.Log.Base" <- function(x, dev)
{
    if (.check.grDevice() && dev) {
        dev.new(title = "Log Likelihood Traceplots")
    }
    par(mfrow = c(3, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    mixlik      <- x@log$mixlik
    plot(mixlik, type = "l", axes = F,
         col = "gray20", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "mixlik", cex = 0.6,
          line = 3)
    mixprior    <- x@log$mixprior
    plot(mixprior, type = "l", axes = F,
         col = "gray47", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "mixprior", cex = 0.6,
          line = 3)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
    cdpost      <- x@log$cdpost
    plot(cdpost, type = "l", axes = F,
         col = "gray47", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "cdpost", cex = 0.6,
          line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Histograms
### Histograms Poisson: Plots the histograms for the Poisson
### parameters and the weights. 
".hist.Poisson.Base" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms")
    }
    lambda <- x@par$lambda
    weight <- x@weight		
    vars <- cbind(lambda, weight[, seq(1, K - 1)])
    lab.names <- vector("list", 2 * K - 1)
    for (k in seq(1, K)) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in seq(K + 1, 2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }  
    .symmetric.Hist(vars, lab.names)
}

### Densities
### Densities Poisson: Plots Kernel densities for the Poisson
### parameters and the weights.
".dens.Poisson.Base" <- function(x, dev)
{
    K   <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities")
    }
    lambda      <- x@par$lambda
    weight      <- x@weight
    vars        <- cbind(lambda, weight[, seq(1, K - 1)])
    lab.names   <- vector("list", 2 * K - 1)
    for (k in seq(1, K)) {
        lab.names[[k]]  <- bquote(lambda[.(k)])
    }
    for (k in seq(K + 1, 2 * K - 1)) {
        lab.names[[k]]  <- bquote(eta[.(k - K)])
    }
    .symmetric.Dens(vars, lab.names)
}

### Subseq: Creates a subsequence of an MCMC sample.
".subseq.Base" <- function(obj, index)
{
    M               <- dim(obj@weight)[1]
    K               <- dim(obj@weight)[2]
    newM            <- sum(index)
    obj@log$cdpost  <- array(obj@log$cdpost[index],
                              dim = c(newM,1))
    obj@weight      <- obj@weight[index, ]
    obj@entropy     <- array(obj@entropy[index],
                              dim = c(newM, 1))
    obj@ST          <- array(obj@ST[index],
                              dim = c(newM, 1)) 
    ## Check which S stay ##
    storeS  <- ifelse(!all(is.na(obj@S)), dim(obj@S)[2], 0)
    if (storeS != 0) {
        ms      <- M - storeS
        index.S <- index[(ms + 1):M]
        N       <- dim(obj@S)[1]
        if (any(index.S)) {
            obj@S       <- array(obj@S[,index.S], dim = c(N, storeS))
        } else {
            obj@S       <- as.array(NA)
        }
    }
    obj@NK          <- obj@NK[index, ]                             
    return(obj)
}

### swapElements: Permutes the elements in an MCMC sample
### for each row.
".swapElements.Base" <- function(obj, index)
{
    ## Rcpp::export 'swap_cc()'
    obj@weight      <- swap_cc(obj@weight, index)
    ## Rcpp::export 'swapInd_cc()'
    M               <- obj@M
    K               <- ncol(index)
    storeS          <- ifelse(!all(is.na(obj@S)), dim(obj@S)[2], 0)
    if (storeS != 0) {
        index.S     <- matrix(index[(M - storeS + 1):M, ], 
                              ncol = K, byrow = TRUE)
        obj@S       <- swapInd_cc(obj@S, index.S)
    }
    ## Rcpp::export 'swapST_cc()'
    obj@ST          <- swapST_cc(obj@ST, index)
    ## Rcpp::export 'swap_cc()'
    obj@NK          <- swapInteger_cc(obj@NK, index)
    return(obj)
}
