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
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

setClass("mcmcoutputfixhier",
    representation(
	    hyper = "list"),
	contains = c("mcmcoutputfix"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

setMethod("show", "mcmcoutputfixhier", 
          function(object) {
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), 
                  "\n")
              cat("     M           :", object@M, "\n")
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

setMethod("plot", signature(x = "mcmcoutputfixhier", 
	y = "ANY"), function(x, y = TRUE, ...) {
	if (x@model@dist == "poisson") {
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
	
		## log ##
		if (.check.grDevice() && y) {
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
})

setMethod("plotHist", signature(x = "mcmcoutputfixhier", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
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
})

## Generic defined int 'mcmcoutputfix.R' ##
setMethod("subseq", signature(object = "mcmcoutputfixhier",
                              index = "array"), 
          function(object, index) {
              ## TODO: Check arguments via .validObject ##
              if (dim(index)[1] != object@M) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "logical") {
                  stop("Argument 'index' must be of type 'logical'.")
              }
              dist <- object@model@dist
              
              ## Call 'subseq()' from 'mcmcoutputfix'
              object <- callNextMethod(object, index)
              
              ## hyper ##
              if (dist == "poisson") {
                  object@hyper$b <- matrix(object@hyper$b[index], 
                                           nrow = object@M, ncol = 1)
              }
              return(object)
          }
)

## Generic defined in 'mcmcoutputfix.R' ##
setMethod("swapElements", signature(object = "mcmcoutputfixhier", 
                                    index = "array"),
          function(object, index) {
              ## Check arguments, TODO: .validObject ##
              if (dim(index)[1] != object@M || dim(index)[2] != object@model@K) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "integer") {
                  stop("Argument 'index' must be of type 'integer'")
              }
              if (!all(index > 0) || any(index > object@model@K)) {
                  stop("Elements in argument 'index' must be greater 0 
                       and must not exceed its number of columns.")
              }
              if (object@model@K == 1) {
                  return(object)
              } else {
                  ## Call method 'swap()' from 'mcmcoutputfix' 
                  object <- callNextMethod(object, index)
                  return(object)
              }
          }
)

setGeneric("getHyper", function(object) standardGeneric("getHyper"))
setMethod("getHyper", "mcmcoutputfixhier", function(object) {
							return(object@hyper)
						}
)
	
## No setters for this object as it is not intended 	##
## that users manipulate this object 		    	##
