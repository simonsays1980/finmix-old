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

.mcmcoutputpermfix <- setClass("mcmcoutputpermfix",
                               contains = c("mcmcpermfix", "mcmcoutputfix"),
                               validity = function(object) 
                               {
                                   ## else: OK
                                   TRUE
                               }
)

setMethod("initialize", "mcmcoutputpermfix", 
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), logperm = list()) 
          {
              .Object@M         <- mcmcoutput@M
              .Object@burnin    <- mcmcoutput@burnin
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object@logperm   <- logperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfix",
          function(object) 
          {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plotTraces", signature(x     = "mcmcoutputpermfix", 
                                  dev   = "ANY",
                                  lik   = "ANY"), 
	function(x, dev = TRUE, lik = 1, ...) 
    {
        if (lik %in% c(0, 1)) {
            if(x@model@dist == "poisson") {
                .permtraces.Poisson(x, dev)
            }
        }
        if (lik %in% c(1, 2)) {
            ## log ##
            .permtraces.Log(x, dev)
        }
	}
)

setMethod("plotHist", signature(x   = "mcmcoutputpermfix", 
                                dev = "ANY"), 
	function(x, dev = TRUE, ...) 
    {
        if(x@model@dist == "poisson") {
            .permhist.Poisson(x, dev)
        }
    }
)

setMethod("plotDens", signature(x   = "mcmcoutputpermfix",
                                dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permdens.Poisson(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputpermfix",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if (x@model@dist == "poisson") {
                  .permpointproc.Poisson(x, dev)
              }
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputpermfix",
                                   dev  = "ANY"),
          function(x, dev, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permsamprep.Poisson(x, dev)
              }
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputpermfix",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permpostdens.Poisson(x, dev)
              }
          }
)

### Private functions.
### These functions are not exported.

### Traces
### Traces Poisson: Plots the traces of the Poisson parameter.
".permtraces.Poisson" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K
    if (.check.grDevice() && dev) {	
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0), 
        oma = c(4, 5, 4,4))
    lambda <- x@parperm$lambda
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

### Traces log-likelihood: Plots the traces of the log-
### likelihoods.
".permtraces.Log" <- function(x, dev)
{
    if(.check.grDevice() && dev) {
        dev.new(title = "Log Likelihood Traceplots")
    }
    par(mfrow = c(2, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    mixlik <- x@logperm$mixlik
    plot(mixlik, type = "l", axes = F,
         col = "gray20", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "mixlik", cex = 0.6,
          line = 3)
    mixprior <- x@logperm$mixprior
    plot(mixprior, type = "l", axes = F,
         col = "gray47", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 3, "mixprior", cex = 0.6,
          line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Histograms
### Histograms Poisson: Plots histograms for all Poisson
### parameters.
".permhist.Poisson" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    lambda <- x@parperm$lambda
    lab.names <- vector("list", K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Hist(lambda, lab.names)
}

### Densities
### Densities Poisson: Plots densities for all Poisson
### parameters.
".permdens.Poisson" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities (permuted)")
    }
    lambda <- x@parperm$lambda
    lab.names <- vector("list", K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Dens(lambda, lab.names)
}

### Plot Point Processes
### Plot Point Process Poisson: Plots the point process
### for the MCMC draws for lambda. The values are plotted
### against a random normal sample. 
".permpointproc.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    M   <- x@M
    if (.check.grDevice() && dev) {
        dev.new("Point Process Representation (MCMC permuted)")
    }
    y.grid  <- replicate(K, rnorm(M))
    if (median(x@parperm$lambda) < 1) {
        lambda  <- log(x@parperm$lambda)
    } else {
        lambda  <- x@parperm$lambda
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
".permsamprep.Poisson" <- function(x, dev)
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
    lambda  <- x@parperm$lambda
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
".permpostdens.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    if (K != 2) {
        warning(paste("A plot of the posterior density is ",
                      "available only for K = 2.", sep = ""))
    } else {
        M   <- x@M
        n   <- min(2000, M)
        lambda  <- x@parperm$lambda
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
