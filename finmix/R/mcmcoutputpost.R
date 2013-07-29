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

setClass("mcmcoutputpost",
	representation(
		post = "list"),
	contains = c("mcmcoutputbase"),
	validity = function(object) {
		## else: OK
		TRUE
	}
)

setMethod("show", "mcmcoutputpost", function(object) {
	cat("Object 'mcmcoutputpost\n")
	cat("	class		:", class(object), "\n")
	cat("	M		:", object@M, "\n")
	cat("	ranperm		:", object@ranperm, "\n")
	cat("	par		: List of", 
		length(object@par), "\n")
	cat("	weight		:", paste(dim(object@weight), 
		collapse = "x"), "\n")
	cat("	log		: List of ",
		length(object@par), "\n")
	cat("	post		: List of ", 
		length(object@post), "\n")
	cat("	entropy:	:", paste(dim(object@entropy),
		collapse = "x"), "\n")
	cat("	ST		:", paste(dim(object@ST),
		collapse = "x"), "\n")
	cat("	S		:", paste(dim(object@S),
		collapse = "x"), "\n")
	cat("	NK		:", paste(dim(object@NK),
		collapse = "x"), "\n")
	cat("	clust		:", paste(dim(object@clust),
		collapse = "x"), "\n")
	cat("	model		: Object of class",
		class(object@model), "\n")
	cat("	prior		: Object of class",
		class(object@prior), "\n")
})

setMethod("show", "mcmcoutputpost", 
          function(object){
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of ", 
                  length(object@par), "\n")
              cat("     weight      : ",
                  paste(dim(object@weight), collapse = "x"), "\n")
              cat("     log         : List of ", 
                  length(object@log), "\n")
              cat("     post        : List of ",
                  length(object@post), "\n")
              cat("     ST          : ", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              cat("     S           : ", 
                  paste(dim(object@S), collapse = "x"), "\n")
              cat("     NK          : ",
                  paste(dim(object@NK), collapse = "x"), "\n")
              cat("     clust       : ",
                  paste(dim(object@clust), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plot", signature(x = "mcmcoutputpost", 
	y = "missing"), function(x, y, ...) {
	if (x@model@dist == "poisson") {
		K <- x@model@K
		trace.n <- K * 2 -1
		if (.check.grDevice()) {
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

setMethod("plotHist", signature(x = "mcmcoutputpost", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K 
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@par$lambda
        weight <- x@weight
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

setMethod("subseq", signature(object = "mcmcoutputpost", 
                              index = "array"), 
          function(object, index) {
              ## TODO: Check arguments via .validObject ##
              if (dim(index)[1] != object@M) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "logical") {
                  stop("Argument 'index' must be of type 'logical'.")
              }
              ## Call 'subseq()' method from 'mcmcoutputfixpost'
              ## class
              object <- callNextMethod(object, index)              
              ## Change owned slots ##
              dist <- object@model@dist
              if (dist == "poisson") {
                  if (object@model@K == 1) {
                      object@post$par$a <- matrix(object@post$par$a[index],
                                                  nrow = object@M, ncol = 1)
                      object@post$par$b <- matrix(object@post$par$b[index],
                                                  nrow = object@M, ncol = 1)
                  } else {
                      object@post$par$a     <- object@post$par$a[index,]
                      object@post$par$b     <- object@post$par$b[index,]
                      object@post$weight    <- object@post$weight[index,]
                  }
              }
              return(object)
          }
)

## Generic defined in 'mcmcoutputfix.R' ##
setMethod("swapElements", signature(object = "mcmcoutputpost", 
                                    index = "array"),
          function(object, index) {
              ## Check arguments, TODO: .validObject ##
              if (dim(index)[1] != object@M || dim(index)[2] != object@model@K) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "integer") {
                  stop("argument 'index# must be of type 'integer'.")
              }
              if (!all(index > 0) || any(index > object@model@K)) {
                  stop("Elements in argument 'index' must be greater 0 
                       and must not exceed its number of columns.")
              }
              if( object@model@K == 1) {
                  return(object) 
              } else {
                  dist <- object@model@dist
                  ## Call method 'swapElements()' from 'mcmcoutputbase' 
                  object <- callNextMethod(object, index)
                  if (dist == "poisson") {
                      ## Rcpp::export 'swap_cc()'
                      object@post$par$a <- swap_cc(object@post$par$a, 
                                                index)
                      object@post$par$b <- swap_cc(object@post$par$b,
                                                index)
                  }
                  return(object)
              }
          }
)

setMethod("getPost", "mcmcoutputpost", function(object) {
						return(object@post)
					}
)
