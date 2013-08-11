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

setClass("mcmcoutputbase", 
         representation(weight 	= "array",
		                entropy	= "array",
                		ST 	    = "array",
                		S 	    = "array",
                		NK 	    = "array",
                		clust 	= "array"),
                    	contains = c("mcmcoutputfix"),
       	 validity = function(object) {
		    ## else: OK
		    TRUE
	     }
)

setMethod("show", "mcmcoutputbase", 
          function(object){
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              cat("     S           :", 
                  paste(dim(object@S), collapse = "x"), "\n")
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
	y = "ANY"), function(x, yi = TRUE, ...) {
	if (x@model@dist == "poisson") {
		K <- x@model@K
		trace.n <- K * 2 - 1
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

setMethod("plotHist", signature(x = "mcmcoutputbase", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K 
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@par$lambda
        weight <- x@weight		
        vars <- cbind(lambda, weight[, seq(1:(K - 1))])
        lab.names <- vector("list", 2 * K - 1)
        for (k in 1:K) {
            lab.names[[k]] <- bquote(lambda[.(k)])
        }
        for (k in (K + 1):(2 * K - 1)) {
            lab.names[[k]] <- bquote(eta[.(k - K)])
        }
        .symmetric.Hist(vars, lab.names)
	}
	
})

## Generic defined in 'mcmcoutputfix.R' ##
setMethod("subseq", signature(object = "mcmcoutputbase", 
                              index = "array"), 
          function(object, index) {
              ## TODO: Check arguments via .validObject ##
              if (dim(index)[1] != object@M) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "logical") {
                  stop("Argument 'index' must be of type 'logical'.")
              }
              M <- object@M
              ## Call 'subseq()' method from 'mcmcoutputfix'
              object <- callNextMethod(object, index)
              ## Change owned slots ##
              object@log$cdpost     <- matrix(object@log$cdpost[index],
                                              nrow = object@M, ncol = 1)
              object@weight         <- object@weight[index, ]
              object@entropy        <- matrix(object@entropy[index],
                                               nrow = object@M, ncol = 1)
              object@ST             <- matrix(object@ST[index], 
                                              nrow = object@M, ncol = 1) 
              ## Check which S stay ##
              stores <- dim(object@S)[2]
              ms <- M - stores
              index.S <- index[(ms + 1):M]
              if(any(index.S) && stores != 0) {
                  object@S <- object@S[,index.S]
              }
              object@NK             <- object@NK[index,]
              return(object)
          }
)

## Generic defined in 'mcmcoutputfix.R' ##
setMethod("swapElements", signature(object = "mcmcoutputbase", 
                                    index = "array"),
          function(object, index) {
              ## Check arguments, TODO: .validObject ##
              if (dim(index)[1] != object@M || dim(index)[2] != object@model@K) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "integer") {
                  stop("Argument 'index' must be of type 'integer'.") 
              }
              if (!all(index > 0) || any(index > object@model@K)) {
                  stop("Elements of argument 'index' must be greater 0 
                       and must not exceed its number of columns.")
              }
              dist          <- object@model@dist
              ## Call method 'swapElements()' from 'mcmcoutputfix' 
              object        <- callNextMethod(object, index)
              ## Rcpp::export 'swap_cc()'
              object@weight <- swap_cc(object@weight, index)
              ## Rcpp::export 'swapInd_cc()'
              M             <- object@M
              storeS        <- dim(object@S)[2]
              if (storeS != 0) {
                  index.S       <- index[(M - storeS + 1):M, ]
                  object@S      <- swapInd_cc(object@S, index.S)
              }
              ## Rcpp::export 'swapST_cc()'
              object@ST     <- swapST_cc(object@ST, index)
              ## Rcpp::export 'swap_cc()'
              object@NK     <- swapInteger_cc(object@NK, index)
              return(object)
          }
)

setGeneric("getWeight", function(object) standardGeneric("getWeight"))
setMethod("getWeight", "mcmcoutputbase", function(object) {
							return(object@weight)
						}
)
setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))
setMethod("getEntropy", "mcmcoutputbase", function(object) {
							return(object@entropy)	
						}
)
setGeneric("getST", function(object) standardGeneric("getST"))
setMethod("getST", "mcmcoutputbase", function(object) {
							return(object@ST)	
						}
)
setGeneric("getS", function(object) standardGeneric("getS"))
setMethod("getS", "mcmcoutputbase", function(object) {
							return(object@S)	
						}
)
## Generic set in 'groupmoments.R' ##
setMethod("getNK", "mcmcoutputbase", function(object) {
							return(object@NK)	
						}
)
setGeneric("getClust", function(object) standardGeneric("getClust"))
setMethod("getClust", "mcmcoutputbase", function(object) {
							return(object@clust)	
						}
)


