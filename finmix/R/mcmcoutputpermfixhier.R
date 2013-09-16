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

.mcmcoutputpermfixhier <- setClass("mcmcoutputpermfixhier",
                                   contains = c("mcmcpermfix", "mcmcoutputfixhier"),
                                   validity = function(object) 
                                   {
                                       ## else: OK
                                       TRUE
                                   }
)

setMethod("initialize", "mcmcoutputpermfixhier",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), logperm = list()) 
          {
              .Object@M         <- mcmcoutput@M
              .Object@burnin    <- mcmcout@burnin
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@hyper     <- mcmcoutput@hyper
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object@logperm   <- logperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixhier",
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
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
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

setMethod("plot", signature(x = "mcmcoutputpermfixhier", 
                            y = "ANY"), 
          function(x, y = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permtraces.Poisson.Hier(x, y)
              }             
              ## log ##
              .permtraces.Log(x, y)
          }
)

setMethod("plotHist", signature(x   = "mcmcoutputpermfixhier", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .permhist.Poisson.Hier(x, dev)
              }
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputpermfixhier", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .permdens.Poisson.Hier(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputpermfixhier",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if (x@model@dist == "poisson") {
                  .permpointproc.Poisson(x, dev)
              }
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputpermfixhier",
                                   dev  = "ANY"),
          function(x, dev, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permsamprep.Poisson(x, dev)
              }
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputpermfixhier",
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
### Traces Poisson: Plots the traces of Poisson parameters 
### and the hyper-parameter 'b'.
".permtraces.Poisson.Hier" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K + 1
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    lambda <- x@parperm$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]),
              cex = 0.6, line = 3)
    }
    b <- x@hyper$b
    plot(b, type = "l", axes = F, 
         col = "gray68", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Histograms
### Histograms Poisson: Plots histograms for all Poisson 
### parameters and the hyper-parameter 'b'.
".permhist.Poisson.Hier." <- function(x, dev)
{
    K <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
	}
    lambda  <- x@parperm$lambda
    b       <- x@hyper$b
    vars    <- cbind(lambda, b)
    lab.names <- vector("list", K + 1)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    lab.names[[K + 1]] <- "b"
    .symmetric.Hist(vars, lab.names)
}

### Densities
### Densities Poisson: Plots densities for all Poisson 
### parameters and the hyper-parameter 'b'.
".permdens.Poisson.Hier." <- function(x, dev)
{
    K <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
	}
    lambda  <- x@parperm$lambda
    b       <- x@hyper$b
    vars    <- cbind(lambda, b)
    lab.names <- vector("list", K + 1)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    lab.names[[K + 1]] <- "b"
    .symmetric.Dens(vars, lab.names)
}
