setClass("mcmcpermfix", 
         representation(
                        Mperm       = "integer",
                        parperm     = "list",
                        logperm     = "list"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

## Define Classes in inheritance structure ##
setClass("mcmcoutputpermfix",
         contains = c("mcmcpermfix", "mcmcoutputfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermfix", 
          function(.Object, mcmcoutput, Mperm, parperm,
                   logperm) {
              .Object@M         <- mcmcoutput@M
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
          function(object) {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
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

setMethod("plot", signature(x = "mcmcoutputpermfix", y = "missing"), 
	function(x, y, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
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
		
		## log ##
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
})

setMethod("plotHist", signature(x = "mcmcoutputpermfix", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K 
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@parperm$lambda
		if (K == 1) {
			hist(lambda, col = "gray65", border = "white",
				cex = 0.7, cex.axis = 0.7, freq = TRUE,
				xlab = "", main = "")
			rug(lambda, col = "gray47")
			mtext(side = 1, bquote(lambda), cex = 0.7,
				line = 3)
		}
		else if (K == 2 || K == 3) {
			par(mfrow = c(1, K), mar = c(2, 2, 2, 2),
				oma = c(4, 5, 1, 5))
			for (k in 1:K) {
				hist(lambda[, k], col = "gray65",
					border = "white", cex = 0.7,
					cex.axis = 0.7, freq = TRUE,
					xlab = "", main = "")
				rug(lambda[, k], col = "gray47")
				mtext(side = 1, bquote(lambda[k = .(k)]),
					cex = 0.7, line = 3)
			}
		}
		else if (K > 3 && K < 17 && sqrt(K)%%1 == 0) {
			par(mfrow = c(sqrt(K), sqrt(K)),
				mar = c(2, 2, 2, 2), 
				oma = c(4, 5, 1, 5))
			for(k in 1:K) {
				hist(lambda[, k], col = "gray65", 
					border = "white", cex = 0.7,
					cex.axis = 0.7, freq = TRUE,
					xlab = "", main = "")
				rug(lambda[, k], col = "gray47")
				mtext(side = 1, bquote(lambda[k = .(k)]),
					cex = 0.7, line = 3)
			}
		}   
		else {
			if(K %% 2 == 0) {
				## check how many rows can be completely 
				## filled
				n <- K %/% 4
				par(mfrow = c(n, 4), mar = c(2, 2, 2, 2),
				oma = c(4, 5, 1, 5))
				for (k in 1:(n * 4)) {
					hist(lambda[, k], col = "gray65", 
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, k], col = "gray47")
					mtext(side = 1, bquote(lambda[k =.(k)]),
						cex = 0.7, line = 3)
				}
				## if some rows cannot be completely filled
				## fill them symmetrically
				if (K %% 4 != 0) {
					## there can be only two left
					hist(lambda[, K - 1], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K - 1], col = "gray47")
					mtext(side = 1, bquote(lambda[k = 
						.(K - 1)]), cex = 0.7, line = 3)
					replicate(2, plot.new())
					hist(lambda[, K], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K], col = "gray47")
					mtext(side = 1, bquote(lambda[k = 
						.(K)]), line = 3)
				}
				
			} else {
				n <- K %/% 5
				par(mfrow = c(n, 5), mar = c(2, 2, 2, 2),
					oma = c(4, 5, 1, 5))
				for (k in 1:(n * 5)) {
					hist(lambda[, k], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, k], col = "gray47")
					mtext(side = 1, bquote(lambda[k = 
						.(k)]), line = 3)
				}
				if (K %% 5 != 1) {
					## put the last one in the middle
					replicate(2, plot.new())
					hist(lambda[, K], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K], col = "gray47")
					mtext(side = 1, bquote(lambda[k = .(K)]),
						line = 3)
					replicate(2, plot.new())		
				}
				else if(K %% 5 != 3) {
					hist(lambda[, K - 2], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K - 2], col = "gray47")
					mtext(side = 1, bquote(lambda[k = 
						.(K - 2)]), line = 3)
					plot.new()
					hist(lambda[, K - 1], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K - 1], col = "gray47")	
					mtext(side = 1, bquote(lambda[k = 
						.(K - 1)]), line = 3)
					plot.new()
					hist(lambda[, K], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K], col = "gray47")
					mtext(side = 1, bquote(lambda[k = 
						.(K)]), line = 3)
				}
			}
		}
	}
	
})

### --- Class 'mcmcoutputpermfixhier' --- ###
setClass("mcmcoutputpermfixhier",
         contains = c("mcmcpermfix", "mcmcoutputfixhier"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermfixhier",
          function(.Object, mcmcoutput, Mperm, parperm, 
                   logperm) {
              .Object@M         <- mcmcoutput@M
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
          function(object) {
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
	y = "ANY"), function(x, y = TRUE, ...) {
	if (x@model@dist == "poisson") {
		K <- x@model@K
		trace.n <- K + 1
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

setMethod("plotHist", signature(x = "mcmcoutputpermfixhier", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@parperm$lambda
        b <- x@hyper$b
	    if (K + 1 < 4) {
			par(mfrow = c(1, K + 1), mar = c(2, 2, 2, 2),
				oma = c(4, 5, 1, 5))
			for (k in 1:K) {
				hist(lambda[, k], col = "gray65",
					border = "white", cex = 0.7,
					cex.axis = 0.7, freq = TRUE,
					xlab = "", main = "")
				rug(lambda[, k], col = "gray47")
				mtext(side = 1, bquote(lambda[k = .(k)]),
					cex = 0.7, line = 3)
			}
            hist(b, col = "gray65", border = "white",
                cex = 0.7, cex.axis = 0.7, freq  = TRUE,
                xlab = "", main ="")
            rug(b, col = "gray47")
            mtext(side = 1, "b", line = 3)
		}
		else if (K + 1> 3 && K + 1 < 17 && sqrt(K + 1)%%1 == 0) {
			par(mfrow = c(sqrt(K + 1), sqrt(K + 1)),
				mar = c(2, 2, 2, 2), 
				oma = c(4, 5, 1, 5))
			for(k in 1:K) {
				hist(lambda[, k], col = "gray65", 
					border = "white", cex = 0.7,
					cex.axis = 0.7, freq = TRUE,
					xlab = "", main = "")
				rug(lambda[, k], col = "gray47")
				mtext(side = 1, bquote(lambda[k = .(k)]),
					cex = 0.7, line = 3)
			}
            hist(b, col = "gray65", border = "white",
                cex = 0.7, cex.axis = 0.7, freq  = TRUE,
                xlab = "", main ="")
            rug(b, col = "gray47")
            mtext(side = 1, "b", line = 3)
		}   
		else {
			if((K + 1) %% 2 == 0) {
				## check how many rows can be completely 
				## filled
				n <- (K + 1) %/% 4
				par(mfrow = c(n, 4), mar = c(2, 2, 2, 2),
				oma = c(4, 5, 1, 5))
				for (k in 1:(n * 4)) {
                    if (k == K + 1) {
                         hist(b, col = "gray65", border = "white",
                            cex = 0.7, cex.axis = 0.7, freq  = TRUE,
                            xlab = "", main ="")
                         rug(b, col = "gray47")
                         mtext(side = 1, "b", line = 3)
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
				if ((K + 1) %% 4 != 0) {
					## there can be only two left
					hist(lambda[, K], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K], col = "gray47")
					mtext(side = 1, bquote(lambda[K = 
						.(K)]), cex = 0.7, line = 3)
					replicate(2, plot.new())
					hist(b, col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(b, col = "gray47")
					mtext(side = 1, "b", line = 3)
				}
				
			} else {
				n <- (K + 1) %/% 5
				par(mfrow = c(n, 5), mar = c(2, 2, 2, 2),
					oma = c(4, 5, 1, 5))
				for (k in 1:(n * 5)) {
                    if (k == (K + 1)) {
                         hist(b, col = "gray65", border = "white",
                            cex = 0.7, cex.axis = 0.7, freq  = TRUE,
                            xlab = "", main ="")
                         rug(b, col = "gray47")
                         mtext(side = 1, "b", line = 3)
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
				if ((K + 1) %% 5 != 1) {
					## put the last one in the middle
					replicate(2, plot.new())
					hist(b, col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(b, col = "gray47")
					mtext(side = 1, "b", line = 3)
					replicate(2, plot.new())		
				}
				else if((K + 1) %% 5 != 3) {
					hist(lambda[, K - 1], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K - 1], col = "gray47")
                    mtext(side = 1, bquote(lambda[k = 
						.(K - 1)]), line = 3)
					plot.new()
					hist(lambda[, K], col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(lambda[, K], col = "gray47")	
					mtext(side = 1, bquote(lambda[k = 
						.(K)]), line = 3)
					plot.new()
					hist(b, col = "gray65",
						border = "white", cex = 0.7,
						cex.axis = 0.7, freq = TRUE,
						xlab = "", main = "")
					rug(b, col = "gray47")
					mtext(side = 1, "b", line = 3)
				}
			}
		}
	}
	
})

## Getters ##
setGeneric("getMperm", function(object) standardGeneric("getMperm"))
setMethod("getMperm", "mcmcpermfix", 
          function(object) {
              w
              return(object@Mperm)
          }
)

setGeneric("getParperm", function(object) standardGeneric("getParperm"))
setMethod("getParperm", "mcmcpermfix", 
          function(object) {
              return(object@parperm)
          }
)

setGeneric("getLogperm", function(object) standardGeneric("getLogperm"))
setMethod("getLogperm", "mcmcpermfix", 
          function(object) {
              return(object@logperm)
          }
)
## No setters as users are not intended to modify these ##
## obect.                                               ##
