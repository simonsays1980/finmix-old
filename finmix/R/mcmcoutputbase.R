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
              cat("     par         : List of ", 
                  length(object@par), "\n")
              cat("     log         : List of ", 
                  length(object@log), "\n")
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
    				rug(lambda, col = "gray47")
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
              stores <- ncol(object@S)
              ms <- M - stores
              index.S <- index[(ms + 1):M]
              if(!any(index.S)) {
                  object@S <- array()
              } else {
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
              storeS        <- ncol(object@S)
              index.S       <- index[(M - storeS + 1):M, ]
              object@S      <- swapInd_cc(object@S, index.S)
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
setGeneric("getNK", function(object) standardGeneric("getNK"))
setMethod("getNK", "mcmcoutputbase", function(object) {
							return(object@NK)	
						}
)
setGeneric("getClust", function(object) standardGeneric("getClust"))
setMethod("getClust", "mcmcoutputbase", function(object) {
							return(object@clust)	
						}
)


