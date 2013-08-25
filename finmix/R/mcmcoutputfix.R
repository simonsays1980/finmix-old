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

.mcmcoutputfix <- setClass("mcmcoutputfix",
                           representation(M 	    = "integer",
                                          ranperm   = "logical",
                                          par 	    = "list",
                                          log	    = "list",
                                          model 	= "model",
                                          prior	    = "prior"),
                           validity = function(object) 
                           {
                               ##else: OK
                               TRUE
                           },
                           prototype(M          = integer(),
                                     ranperm    = logical(),
                                     par        = list(),
                                     log        = list(),
                                     model      = model(),
                                     prior      = prior()
                                     )
)

setMethod("show", "mcmcoutputfix", 
          function(object) {
              cat("Object 'mcmcoutputfix'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plot", signature(x = "mcmcoutputfix", 
                            y = "missing"), 
          function(x, y, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .traces.Poisson(x, dev)
                  ## log ##
                  .traces.Log(x, dev)
              }
          }
)

setMethod("plotHist", signature(x   = "mcmcoutputfix", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .hist.Poisson(x, dev)
              }
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputfix",
                                dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if(x@model@dist == "poisson") {
                  .dens.Poisson(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputfix",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if (x@model@dist == "poisson") {
                  .pointproc.Poisson(x, dev)
              }
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputfix",
                                   dev  = "ANY"),
          function(x, dev, ...) 
          {
              if (x@model@dist == "poisson") {
                  .samprep.Poisson(x, dev)
              }
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputfix",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .postdens.Poisson(x, dev)
              }
          }
)

setMethod("subseq", signature(object    = "mcmcoutputfix", 
                              index     = "array"), 
          function(object, index) 
          {
              .subseq.valid.Arg(object, index)
              dist      <- object@model@dist
              object@M  <- sum(index)
              ## log ##
              object    <- .subseq.Log.Fix(object, index)
              ## par ##
              if(dist == "poisson") {
                  .subseq.Poisson(object, index)
              }
          }
)

setMethod("swapElements", signature(object  = "mcmcoutputfix", 
                                    index   = "array"),
          function(object, index) 
          { ## Check arguments, TODO: .validObject ##
              .swapElements.valid.Arg(object, index)
              if(object@model@K == 1) {
                  return(object)
              } else {
                  dist <- object@model@dist
                  if (dist == "poisson") {
                      .swapElements.Poisson(object, index)
                  }
              }
          }
) 
              
## Getters ##
setMethod("getM", "mcmcoutputfix", 
          function(object) 
          {
              return(object@M)
          }
)

setMethod("getRanperm", "mcmcoutputfix", 
          function(object) 
          {
              return(object@ranperm)
          }
)

setMethod("getPar", "mcmcoutputfix", 
          function(object) 
          {
              return(object@par)
          }
)

setMethod("getLog", "mcmcoutputfix", 
          function(object) 
          {
              return(object@log)
          }
)

setMethod("getModel", "mcmcoutputfix", 
          function(object) 
          {
              return(object@model)
          }
)

setMethod("getPrior", "mcmcoutputfix", 
          function(object) 
          {
              return(object@prior)
          }
)

## No setters as users are not intended to manipulate ##
## this object ##i

### Private functions
### These functions are not exported

### Plot
### Traces Poisson: Plots the traces of MCMC samples
### for Poisson mixture. If dev = FALSE, no graphical
### device is started, instead it is assumed that the
### user wants to save the graphic to a file.
".traces.Poisson" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K
    if (.check.grDevice() && dev) {	
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0), 
        oma = c(4, 5, 4,4))
    lambda <- x@par$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]), 
              cex = 0.6, line = 3)
    }
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)		
}

### Traces Poisson: Plots the traces of MCMC samples
### for the log-likelihoods. If dev = FALSE, no graphical
### device is started, instead it is assumed that the
### user wants to save the graphic to a file.
".traces.Log" <- function(x, dev)
{
    if(.check.grDevice() && dev) {
        dev.new(title = "Log Likelihood Traceplots")
    }
    par(mfrow = c(2, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    mixlik <- x@log$mixlik
    plot(mixlik, type = "l", axes = F,
         col = "gray20", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "mixlik", cex = 0.6,
          line = 3)
    mixprior <- x@log$mixprior
    plot(mixprior, type = "l", axes = F,
         col = "gray47", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "mixprior", cex = 0.6,
          line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)		
}

### Plot Histogramms
### Plot Hist Poisson: Plots Histograms for each component
### parameter. 
".hist.Poisson" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms")
    }
    lambda <- x@par$lambda
    if (K == 1) {
        .symmetric.Hist(lambda, list(bquote(lambda)))
    } else {
        lab.names <- vector("list", K)
        for (k in 1:K) {
            lab.names[[k]] <- bquote(lambda[.(k)])
        }
        .symmetric.Hist(lambda, lab.names)
    }
}

### Plot Densities
### Plot Dens Poisson: Plots Kernel densities for each
### component parameter.
".dens.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities")
    }
    lambda  <- x@par$lambda
    if (K == 1) {
        .symmetric.Dens(lambda, list(bquote(lambda)))
    } else {
        lab.names   <- vector("list", K)
        for (k in seq(1, K)) {
            lab.names[[k]]  <- bquote(lambda[.(k)])
        }
        .symmetric.Dens(lambda, lab.names)
    }
}

### Plot Point Processes
### Plot Point Process Poisson: Plots the point process
### for the MCMC draws for lambda. The values are plotted
### against a random normal sample. 
".pointproc.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    M   <- x@M
    if (.check.grDevice() && dev) {
        dev.new("Point Process Representation (MCMC)")
    }
    y.grid  <- replicate(K, rnorm(M))
    if (median(x@par$lambda) < 1) {
        lambda  <- log(x@par$lambda)
    } else {
        lambda  <- x@par$lambda
    }
    col.grid <- gray.colors(K, start = 0, 
                            end = 0.5)
    legend.names    <- vector("list", K)
    for (k in seq(1, K)) {
        legend.names[[k]]   <- bquote(lambda[.(k)])
    }
    plot(lambda, y.grid, pch = 20, col = col.grid,
         cex = .7, cex.axis = .7, cex.lab = .7,
         main = "", ylab = "", xlab = "")
    mtext(side = 1, bquote(lambda), cex = .7, 
          cex.lab = .7, line = 3)
    legend("topright", legend = do.call(expression, 
                                        legend.names),
           col = col.grid, fill = col.grid)
}

### Plot sampling representation
### Plot sampling representation Poisson: Plots the sampling
### representation for Poisson parameters. Each parameter sample
### is combined with the other samples. 
".samprep.Poisson" <- function(x, dev)
{
    K       <- x@model@K
    if (K == 1) {
        warning(paste("Sampling representation is only ",
                      "available for mixture models with ",
                      "K > 1.", sep = ""))
        return(FALSE)
    }
    M       <- x@M
    n       <- min(2000, x@M)
    n.perm  <- choose(K, 2) * factorial(2)
    lambda  <- x@par$lambda
    if (.check.grDevice() && dev) {
        dev.new(title = "Sampling Representation")
    }
    comb    <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
    comb    <- comb[which(comb[, 1] != comb[, 2]), ]
    lambda  <- lambda[seq(1, n), ]
    lambda  <- matrix(lambda[,comb], nrow = n * n.perm, ncol = 2)
    plot(lambda, col = "gray47", cex.lab = .7, cex.axis = .7,
         cex = .7, pch = 20, main = "", xlab = "", ylab = "")
    abline(0, 1, lty = 1)
    mtext(side = 1, bquote(lambda), cex = .7, cex.lab = .7,
          line = 3)
    mtext(side = 2, bquote(lambda), cex = .7, cex.lab = .7,
          line = 3)

}

### Posterior Density
### Posterior Density Poisson: Plots a contour plot of the 
### posterior density of the sampled parameters for K = 2.
".postdens.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    if (K != 2) {
        warning(paste("A plot of the posterior density is ",
                      "available only for K = 2.", sep = ""))
    } else {
        M   <- x@M
        n   <- min(2000, M)
        lambda  <- x@par$lambda
        lambda  <- lambda[seq(1, n), ]
        dens    <- bkde2D(lambda, bandwidth = c(sd(lambda[, 1]),
                                                sd(lambda[, 2])))
        if (.check.grDevice() && dev) {
            dev.new(title = "Posterior Density Contour Plot (MCMC)")
        } 
        contour(dens$x1, dens$x2, dens$fhat, cex = .7, 
                cex.lab = .7, cex.axis = .7, col = "gray47", 
                main = "", xlab = "", ylab = "")
        mtext(side = 1, bquote(lambda[1]), cex = .7, 
              cex.lab = .7, line = 3)
        mtext(side = 2, bquote(lambda[2]), cex = .7,
              cex.lab = .7, line = 3)
        if (.check.grDevice() && dev) {
            dev.new(title = "Posterior Density Persepctive Plot (MCMC)")
        }
        persp(dens$x1, dens$x2, dens$fhat, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, 
              expand = .5, lphi = 180, ltheta = 90, 
              r = 40, d = .1, ticktype = "detailed", zlab = 
              "Density", xlab = "k = 1" , ylab = "k = 2")
    }
}

### Logic
### Logic subseq: This function is used for each 
### distribution type in 'model'. It crreates a subsequence
### for the log-likelihoods. 
".subseq.Log.Fix" <- function(obj, index)
{
    obj@log$mixlik     <- matrix(obj@log$mixlik[index],
                                 nrow = obj@M, ncol = 1)
    obj@log$mixprior   <- matrix(obj@log$mixprior[index],
                                 nrow = obj@M, ncol = 1)
    return(obj)
}

### Logic subseq Poisson: This function creates a subsequence 
### MCMC Poisson parameter samples. 
".subseq.Poisson" <- function(obj, index) 
{
    if(obj@model@K == 1) {
        obj@par$lambda <- matrix(obj@par$lambda[index], 
                                    nrow = obj@M, ncol = 1)
    } else {
        obj@par$lambda <- obj@par$lambda[index,]
    }
    return(obj)
}

### Log swapElements
### Logic swapElements Poisson: This function permutes
### the elements in the MCMC sample for Poisson 
### parameters by calling the C++-function 'swap_cc()'.
".swapElements.Poisson" <- function(obj, index)
{
    ## Rcpp::export 'swap_cc'
    obj@par$lambda <- swap_cc(obj@par$lambda, index)
    return(obj)
}

### Validity
### Validity subseq: The index given to 'subseq()' must 
### have dimension M x 1 and must contain logical values.
".subseq.valid.Arg" <- function(obj, index) 
{
    if (dim(index)[1] != obj@M) {
        stop("Argument 'index' has wrong dimension.")
    }
    if (typeof(index) != "logical") {
        stop("Argument 'index' must be of type 'logical'.")
    }
}

### Validity swapElements: The index given to 'swapElements()'
### must have dimension M x K. It must be of type 'integer'
### and must be in the range 1, ..., K.
".swapElements.valid.Arg" <- function(obj, index)
{
    if (dim(index)[1] != obj@M || dim(index)[2] != obj@model@K) {
        stop("Argument 'index' has wrong dimension.")
    }
    if (typeof(index) != "integer") {
        stop("Argument 'index' must be of type 'integer'.")
    }
    if (!all(index > 0) || any(index > obj@model@K)) {
        stop(paste("Elements in argument 'index' must be greater 0", 
             "and must not exceed its number of columns.", 
             sep = ""))
    }
}
