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

setClass("mcmcoutputfixhierpost", 
         representation(post = "list"),
	     contains = c("mcmcoutputfixhier"),
	     validity = function(object) {
			 ## else: OK
			 TRUE
	     }
)

setMethod("show", "mcmcoutputfixhierpost", function(object) {
	cat("Object 'mcmcoutputfixhierpost\n")
	cat("	class		:", class(object), "\n")
	cat("	M		:", object@M, "\n")
	cat("	ranperm		:", object@ranperm, "\n")
	cat("	par		: List of ", 
		length(object@par), "\n")
	cat("	log		: List of ",
		length(object@par), "\n")
	cat("	hyper		: List of ",
		length(object@hyper), "\n")
	cat("	post		: List of ",
		length(object@post), "\n")
	cat("	model		: Object of class",
		class(object@model), "\n")
	cat("	prior		: Object of class",
		class(object@prior), "\n")
})

setMethod("plot", signature(x = "mcmcoutputfixhierpost", 
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

setMethod("plotHist", signature(x = "mcmcoutputfixhierpost", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@par$lambda
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

setMethod("subseq", signature(object = "mcmcoutputfixhierpost",
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
              ### Call only the 'subseq' method from the 
              ## 'mcmcoutputfixhier' class
              object <- callNextMethod(object, index)
              ## post ##
              if (dist == "poisson") {
                  if (object@model@K == 1) {
                      object@post$par$a <- matrix(object@post$par$a[index], 
                                                  nrow = object@M, ncol = 1)
                      object@post$par$b <- matrix(object@post$par$b[index],
                                                  nrow = object@M, ncol = 1)
                  } else {
                      object@post$par$a <- object@post$par$a[index, ]
                      object@post$par$b <- object@post$par$b[index, ]
                  }
              }
              return(object)
          }
)

## Generic defined in 'mcmcoutputfix.R' ##
setMethod("swapElements", signature(object = "mcmcoutputfixhierpost", 
                                    index = "array"),
          function(object, index) {
              ## Check arguments, TODO: .validObject ##
              if (dim(index)[1] != object@M || dim(index)[2] != object@model@K) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "integer") {
                  stop("Argzment 'index' must be of type 'integer'.")
              }
              if (!all(index > 0) || any(index > object@model@K)) {
                  stop("Elements of argument 'index' must be greater 0 
                       and must not exceed its number of columns.")
              }
              if (object@model@K == 1) {
                  return(object)
              } else {
                  dist <- object@model@dist
                  ## Call method 'swapElements()' from 'mcmcoutputfixhier' 
                  object <- callNextMethod(object, index)
    
                  if (dist == "poisson") {
                      ## Rcpp::export 'swap_cc' 
                      object@post$par$a <- swap_cc(object@post$par$a, index)
                      object@post$par$b <- swap_cc(object@post$par$b, index)
                  }
                  return(object)
              }
          }
)


