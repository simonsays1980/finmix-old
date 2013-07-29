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

setClass("mcmcpermind", 
         representation(weightperm  = "array",
                        entropyperm = "array",
                        STperm      = "array",
                        Sperm       = "array",
                        NKperm      = "array"
                        ),
         contains = c("mcmcpermfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

## Define Classes in inheritance structure ##
setClass("mcmcoutputpermbase",
         contains = c("mcmcpermind", "mcmcoutputbase"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermbase",
          function(.Object, mcmcoutput, Mperm, parperm,  
                   weightperm, logperm, entropyperm,
                   STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@weightperm    <- weightperm
              .Object@logperm       <- logperm
              .Object@entropyperm   <- entropyperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermbase", 
          function(object){
              cat("Object 'mcmcoutputperm'\n")
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
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     weightperm  :",
                  paste(dim(object@weightperm), collapse = "x"), "\n")
              cat("     logperm     : List of", 
                  length(object@logperm), "\n")
              cat("     entropyperm :",
                  paste(dim(object@entropyperm), collapse = "x"), "\n")
              cat("     STperm      :",
                  paste(dim(object@STperm), collapse = "x"), "\n")
              cat("     Sperm       :",
                  paste(dim(object@Sperm), collapse = "x"), "\n")
              cat("     NKperm      :", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plot", signature(x = "mcmcoutputpermbase", 
	y = "ANY"), function(x, yi = TRUE, ...) {
	if (x@model@dist == "poisson") {
		K <- x@model@K
		trace.n <- K * 2 - 1
		if (.check.grDevice() && y) {
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
		weight <- x@weightperm
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

setMethod("plotHist", signature(x = "mcmcoutputpermbase", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K 
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@parperm$lambda
        weight <- x@weightperm
        hist.n <- K * 2 - 1
		if (hist.n == 1) {
			hist(lambda, col = "gray65", border = "white",
				cex = 0.7, cex.axis = 0.7, freq = TRUE,
				xlab = "", main = "")
			rug(lambda, col = "gray47")
			mtext(side = 1, bquote(lambda), cex = 0.7,
				line = 3)
		} else if (hist.n == 3) {
			par(mfrow = c(1, hist.n), mar = c(2, 2, 2, 2),
				oma = c(4, 5, 1, 5))
			for (k in 1:hist.n) {
                if (k > K) {
                    i <- k - K
                    hist(weight[, i], col = "gray65", 
                         border = "white", cex = 0.7,
                         cex.axis = 0.7, freq = TRUE,
                         xlab = "", main = "")
                    rug(weight[, i], col = "gray47")
                    mtext(side = 1, bquote(eta[i = .(i)]),
                          line = 3)
                } else {
    				hist(lambda[, k], col = "gray65",
    					border = "white", cex = 0.7,
    					cex.axis = 0.7, freq = TRUE,
    					xlab = "", main = "")
    				rug(lambda[, k], col = "gray47")
    				mtext(side = 1, bquote(lambda[k = .(k)]),
    					cex = 0.7, line = 3)
                }
			}
		}
		else if (hist.n > 3 && hist.n < 17 && sqrt(hist.n)%%1 == 0) {
			par(mfrow = c(sqrt(hist.n), sqrt(hist.n)),
				mar = c(2, 2, 2, 2), 
				oma = c(4, 5, 1, 5))
			for(k in 1:hist.n) {
                if (k > K) {
                    i <- k - K
                    hist(weight[, i], col = "gray65", 
                        border = "white", cex = 0.7,
                        cex.axis = 0.7, freq = TRUE,
                        xlab = "", main = "")
                    rug(weight[, i], col = "gray47")
                    mtext(side = 1, bquote(eta[i = .(i)]),
                        line = 3)
                } else {
    				hist(lambda[, k], col = "gray65", 
    					border = "white", cex = 0.7,
    					cex.axis = 0.7, freq = TRUE,
    					xlab = "", main = "")
    				rug(lambda[, k], col = "gray47")
    				mtext(side = 1, bquote(lambda[k = .(k)]),
    					cex = 0.7, line = 3)
                }
			}
		}   
		else {
			if(hist.n %% 2 == 0) {
				## check how many rows can be completely 
				## filled
				n <- hist.n %/% 4
				par(mfrow = c(n, 4), mar = c(2, 2, 2, 2),
				oma = c(4, 5, 1, 5))
				for (k in 1:(n * 4)) {
                    if (k > K) {
                        i <- k - K
                        hist(weight[, i], col = "gray65", 
                             border = "white", cex = 0.7,
                             cex.axis = 0.7, freq = TRUE,
                             xlab = "", main = "")
                        rug(weight[, i], col = "gray47")
                        mtext(side = 1, bquote(eta[i = .(i)]),
                            line = 3)
                    } else {
					    hist(lambda[, k], col = "gray65", 
					    	border = "white", cex = 0.7,
					    	cex.axis = 0.7, freq = TRUE,
					    	xlab = "", main = "")
			    		rug(lambda[, k], col = "gray47")
			    		mtext(side = 1, bquote(lambda[k =.(k)]),
			    			cex = 0.7, line = 3)
                    }
				}
				## if some rows cannot be completely filled
				## fill them symmetrically
				if (hist.n %% 4 != 0) {
					## there can be only two left
					hist(weight[, K - 2], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(weight[, K - 2], col = "gray47")
					mtext(side = 1, bquote(eta[k = 
						.(K - 2)]), cex = 0.7, line = 3)
					replicate(2, plot.new())
					hist(weight[, K - 1], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(weight[, K - 1], col = "gray47")
					mtext(side = 1, bquote(eta[k = 
						.(K - 1)]), line = 3)
				}
				
			} else {
				n <- hist.n %/% 5
				par(mfrow = c(n, 5), mar = c(2, 2, 2, 2),
					oma = c(4, 5, 1, 5))
				for (k in 1:(n * 5)) {
                    if (k > K) {
                        i <- k - K
                        hist(weight[, i], col = "gray65", 
                             border = "white", cex = 0.7,
                             cex.axis = 0.7, freq = TRUE,
                             xlab = "", main = "")
                        rug(weight[, i], col = "gray47")
                        mtext(side = 1, bquote(eta[i = .(i)]),
                             line = 3)
                    } else {
    					hist(lambda[, k], col = "gray65",
    						border = "white", cex = 0.7,
    						cex.axis = 0.7, freq = TRUE,
    						xlab = "", main = "")
    					rug(lambda[, k], col = "gray47")
    					mtext(side = 1, bquote(lambda[k = 
    						.(k)]), line = 3)
                    }
				}
				if (hist.n %% 5 != 1) {
					## put the last one in the middle
					replicate(2, plot.new())
					hist(weight[, K - 1], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(weight[, K - 1], col = "gray47")
					mtext(side = 1, bquote(eta[k = .(K - 1)]),
						line = 3)
					replicate(2, plot.new())		
				}
				else if(K %% 5 != 3) {
					hist(weight[, K - 3], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(weight[, K - 3], col = "gray47")
					mtext(side = 1, bquote(eta[k = 
						.(K - 3)]), line = 3)
					plot.new()
					hist(weight[, K - 2], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(weight[, K - 2], col = "gray47")	
					mtext(side = 1, bquote(eta[k = 
						.(K - 2)]), line = 3)
					plot.new()
					hist(weight[, K - 1], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(weight[, K - 1], col = "gray47")
					mtext(side = 1, bquote(eta[k = 
						.(K - 1)]), line = 3)
				}
			}
		}
	}
	
})

### --- Class 'mcmcoutputpermhier' --- ###
setClass("mcmcoutputpermhier",
         contains = c("mcmcpermind", "mcmcoutputhier"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermhier",
          function(.Object, mcmcoutput, Mperm, parperm,
                   weightperm, logperm, entropyperm, 
                   STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@hyper         <- mcmcoutput@hyper
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@weightperm    <- weightperm
              .Object@logperm       <- logperm
              .Object@entropyperm   <- entropyperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermhier", 
          function(object){
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              cat("     S           :", 
                  paste(dim(object@S), collapse = "x"), "\n")
              cat("     NK          :",
                  paste(dim(object@NK), collapse = "x"), "\n")
              cat("     clust       :",
                  paste(dim(object@clust), collapse = "x"), "\n")
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     weightperm  :",
                  paste(dim(object@weightperm), collapse = "x"), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     entropyperm :",
                  paste(dim(object@entropyperm), collapse = "x"), "\n")
              cat("     STperm      :",
                  paste(dim(object@STperm), collapse = "x"), "\n")
              cat("     Sperm       :",
                  paste(dim(object@Sperm), collapse = "x"), "\n")
              cat("     NKperm      :", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plot", signature(x = "mcmcoutputpermhier", 
	y = "missing"), function(x, y, ...) {
	if (x@model@dist == "poisson") {
		K <- x@model@K
		trace.n <- K * 2
		if (.check.grDevice()) {
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
		weight <- x@weightperm
		for (k in 1:(K - 1)) {
			plot(weight[, k], type = "l", axes = F, 
				col = "gray47", xlab = "", ylab = "")
			axis(2, las = 2, cex.axis = 0.7)
			mtext(side = 2, las = 2, bquote(eta[k = .(k)]),
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
		if(.check.grDevice()) {
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

setMethod("plotHist", signature(x = "mcmcoutputpermhier", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K 
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@parperm$lambda
        weight <- x@weightperm
        b <- x@hyper$b
        hist.n <- K * 2
		if (hist.n == 2) {
            par(mfrow = c(1, hist.n), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
			hist(lambda, col = "gray65", border = "white",
				cex = 0.7, cex.axis = 0.7, freq = TRUE,
				xlab = "", main = "")
			rug(lambda, col = "gray47")
			mtext(side = 1, bquote(lambda), cex = 0.7,
				line = 3)
            hist(b, col = "gray65", border = "white",
                cex = 0.7, cex.axis = 0.7, freq = TRUE,
                xlab = "", main = "")
            rug(b, col = "gray47")
            mtext(side = 1, "b", line = 3)
		} else if (hist.n > 3 && hist.n < 17 && sqrt(hist.n) 
                  %% 1 == 0) {
			par(mfrow = c(sqrt(hist.n), sqrt(hist.n)),
				mar = c(2, 2, 2, 2), 
				oma = c(4, 5, 1, 5))
			for(k in 1:hist.n) {
                if (k == hist.n) {
                    hist(b, col = "gray65", border = "white",
                        cex = 0.7, cex.axis = 0.7, freq = TRUE,
                        xlab = "", main = "")
                    rug(b, col = "gray47")
                    mtext(side = 1, "b", line = 3)
                } else if (k > K && k != hist.n) {
                    i <- k - K
                    hist(weight[, i], col = "gray65", 
                        border = "white", cex = 0.7,
                        cex.axis = 0.7, freq = TRUE,
                        xlab = "", main = "")
                    rug(weight[, i], col = "gray47")
                    mtext(side = 1, bquote(eta[i = .(i)]),
                        line = 3)
                } else {
    				hist(lambda[, k], col = "gray65", 
    					border = "white", cex = 0.7,
    					cex.axis = 0.7, freq = TRUE,
    					xlab = "", main = "")
    				rug(lambda[, k], col = "gray47")
    				mtext(side = 1, bquote(lambda[k = .(k)]),
    					cex = 0.7, line = 3)
                }
			}
		}   
		else {
			## check how many rows can be completely 
			## filled
			n <- hist.n %/% 4
			par(mfrow = c(n, 4), mar = c(2, 2, 2, 2),
			oma = c(4, 5, 1, 5))
			for (k in 1:(n * 4)) {
                if (k == hist.n) {
                    hist(b, col = "gray65", border = "white",
                        cex = 0.7, cex.axis = 0.7, freq = TRUE,
                        xlab = "", main = "")
                    rug(b, col = "gray47")
                    mtext(side = 1, "b", line = 3)
                } else if (k > K && k != hist.n) {
                    i <- k - K
                    hist(weight[, i], col = "gray65", 
                        border = "white", cex = 0.7,
                        cex.axis = 0.7, freq = TRUE, 
                        xlab = "", main = "")
                    rug(weight[, i], col = "gray47")
                    mtext(side = 1, bquote(eta[i = .(i)]),
                        line = 3)
                } else {
				    hist(lambda[, k], col = "gray65", 
					    border = "white", cex = 0.7,
					    cex.axis = 0.7, freq = TRUE,
					    xlab = "", main = "")
			    	rug(lambda[, k], col = "gray47")
			    	mtext(side = 1, bquote(lambda[k =.(k)]),
			    		cex = 0.7, line = 3)
                }
			}
			## if some rows cannot be completely filled
			## fill them symmetrically
			if (hist.n %% 4 != 0) {
				## there can be only two left
				hist(weight[, K - 1], col = "gray65",
					border = "white", cex = 0.7,
					cex.axis = 0.7, freq = TRUE,
					xlab = "", main = "")
				rug(weight[, K - 1], col = "gray47")
				mtext(side = 1, bquote(eta[k = 
					.(K - 1)]), cex = 0.7, line = 3)
				replicate(2, plot.new())
				hist(b, col = "gray65",
					border = "white", cex = 0.7,
					cex.axis = 0.7, freq = TRUE,
					xlab = "", main = "")
				rug(b, col = "gray47")
				mtext(side = 1, "b", line = 3)
			}
				
 		}
	}
	
})

## Getters ##
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))
setMethod("getWeightperm", "mcmcpermind", function(object) {
          return(object@weightperm)
})

setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))
setMethod("getEntropyperm", "mcmcpermind", 
          function(object) {
              return(object@entropyperm)
          }
)

setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))
setMethod("getSTperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

setGeneric("getSperm", function(object) standardGeneric("getSperm"))
setMethod("getSperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))
setMethod("getNKperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
