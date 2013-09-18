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

.mcmcoutputfixhier <- setClass("mcmcoutputfixhier",
                            representation(hyper = "list"),
                            contains = c("mcmcoutputfix"),
                            validity = function(object) {
                                ## else: OK
                                TRUE
                            },
                            prototype(hyper = list())
)

setMethod("show", "mcmcoutputfixhier", 
          function(object) 
          {
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), 
                  "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, 
                  "\n")
              cat("     par         : List of",
                  length(object@par), "\n")
              cat("     log         : List of",
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     model       : Object of class",
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plotTraces", signature(x     = "mcmcoutputfixhier", 	
                                  dev   = "ANY",
                                  lik   = "ANY"), 
          function(x, dev = TRUE, lik = 1, ...) 
          {
              if (lik %in% c(0, 1)) {
                  if (x@model@dist == "poisson") {
                      .traces.Poisson.Hier(x, dev)
                  }
              }
              if (lik %in% c(1, 2)) {
                  ## log ##
                  .traces.Log(x, dev)
              }	
          }
)

setMethod("plotHist", signature(x = "mcmcoutputfixhier", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .hist.Poisson.Hier(x, dev)
              }
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputfixhier",
                                dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if (x@model@dist == "poisson") {
                  .dens.Poisson.Hier(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputfixhier",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              ## Call 'plotPointProc()' from 'mcmcoutputfix'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputfixhier",
                                   dev  = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              ## Call 'plotSampRep()' from 'mcmcoutputfix'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputfixhier",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotPostDens()' from 'mcmcoutputfix'
              callNextMethod(x, dev, ...)
          }
)

setMethod("subseq", signature(object = "mcmcoutputfixhier",
                              index = "array"), 
          function(object, index) 
          {
              ## Call 'subseq()' from 'mcmcoutputfix'
              callNextMethod(object, index)
              dist <- object@model@dist             
              ## hyper ##
              if (dist == "poisson") {
                  .subseq.Poisson.Hier(object, index)
              }
          }
)

setMethod("swapElements", signature(object = "mcmcoutputfixhier", 
                                    index = "array"),
          function(object, index) 
          {
              ## Check arguments, TODO: .validObject ##
              .swapElements.valid.Arg(object, index)
              if (object@model@K == 1) {
                  return(object)
              } else {
                  ## Call method 'swap()' from 'mcmcoutputfix' 
                  callNextMethod(object, index)
              }
          }
)

setMethod("getHyper", "mcmcoutputfixhier", 
          function(object) 
          {
              return(object@hyper)
          }
)
	
## No setters for this object as it is not intended ##
## that users manipulate this object. 		    	##

### Private functions
### These functions are not exported.

### Plot

### Plot Traces
### Plot traces Poisson: Plots traces for each component
### parameter of a Poisson mixture and the hyper parameter 'b'.
".traces.Poisson.Hier" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K + 1
    if (.check.grDevice() && y) {
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
    b <- x@hyper$b
    plot(b, type = "l", axes = F, 
         col = "gray68", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)    
    mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)	
}

### Plot Histograms
### Plot hist Poisson: Plots histograms for each component
### parameter and the hyper parameter 'b'.
".hist.Poisson.Hier" <- function(x, dev)
{
    K <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms")
    }
    lambda  <- x@par$lambda
    b       <- x@hyper$b
    vars    <- cbind(lambda, b)
    if (K == 1) {
        lab.names <- list(bquote(lambda), "b")
        .symmetric.Hist(vars, lab.names)
    } else {
        lab.names <- vector("list", K + 1)
        for (k in 1:K) {
            lab.names[[k]] <- bquote(lambda[.(k)])
        }
        lab.names[[K + 1]] <- "b"
        .symmetric.Hist(vars, lab.names)
    }
}

### Plot Densities
### Plot Dens Poisson Hier: Plots Kernel densities for each
### component parameter and the hyper parameter 'b'.
".dens.Poisson.Hier" <- function(x, dev)
{
    K   <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new("Densities")
    }
    lambda  <- x@par$lambda
    b       <- x@hyper$b
    vars    <- cbind(lambda, b)
    if (K == 1) {
        lab.names   <- list(bquote(lambda), "b")
        .symmetric.Dens(vars, lab.names)
    } else {
        lab.names   <- vector("list", K + 1)
        for (k in seq(1, K)) {
            lab.names[[k]]  <- bquote(lambda[.(k)])
        }
        lab.names[[K + 1]]  <- "b"
        .symmetric.Dens(vars, lab.names)
    }
}

### Logic
### Logic subseq Hier: Creates a subsequence for the sample
### of the Poisson hyper parameter 'b'.
".subseq.Poisson.Hier" <- function(obj, index)
{
    obj@hyper$b <- array(obj@hyper$b[index], 
                          dim = c(obj@M, 1))
    return(obj)
}

