setClass("mcmcoutputhierpost", 
         representation(post = "list"),
      	 contains = c("mcmcoutputhier"),
	     validity = function(object) {
			 ## else: OK
			 TRUE
	     }
)

## Set 'mcmcoutput' to the virtual class inheriting 	##
## to each other 'mcmcoutput' class. 			##
## This is done to simplify dispatching methods.	##
setClassUnion("mcmcoutput", 
	c(
		"mcmcoutputfix",
		"mcmcoutputfixhier",
		"mcmcoutputfixpost",
		"mcmcoutputfixhierpost",
		"mcmcoutputbase",
		"mcmcoutputhier",
		"mcmcoutputpost",
		"mcmcoutputhierpost")
)

setMethod("show", "mcmcoutputhierpost", 
          function(object){
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     weight      :",
                  paste(dim(object@weight), collapse = "x"), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     post        : List of",
                  length(object@post), "\n")
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

setMethod("plot", signature(x = "mcmcoutputhierpost", 
	y = "missing"), function(x, y, ...) {
	if (x@model@dist == "poisson") {
		K <- x@model@K
		trace.n <- K * 2 
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

setMethod("plotHist", signature(x = "mcmcoutputhierpost", dev = "ANY"), 
	function(x, dev = TRUE, ...) {
	if(x@model@dist == "poisson") {
		K <- x@model@K 
		if (.check.grDevice() && dev) {
			dev.new(title = "Histograms")
		}
		lambda <- x@par$lambda
        weight <- x@weight
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

setMethod("subseq", signature(object = "mcmcoutputhierpost", 
                              index = "array"), 
          function(object, index) {
              ## TODO: Check arguments via .validObject ##
              if (dim(index)[1] != object@M) {
                  stop("Argument 'index' has wrong dimension.")
              }
              if (typeof(index) != "logical") {
                  stop("Argument 'index' must be of type 'logical'.")
              }
              ## Call 'subseq()' method from 'mcmcoutputfixhierpost'
              ## class
              object    <- callNextMethod(object, index)
              # Change owned slots #
              dist      <- object@model@dist
              K         <- object@model@K
              if (dist == "poisson") {
                  if (K == 1) {
                      object@post$par$a     <- matrix(object@post$par$a[index],
                                                      nrow = object@M, ncol = 1)
                      object@post$par$b     <- matrix(object@post$par$a[index],
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
setMethod("swapElements", signature(object = "mcmcoutputhierpost", 
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
              if (object@model@K == 1) {
                  return(object)
              }
              dist <- object@model@dist
              ## Call method 'swapElements()' from 'mcmcoutputhier' 
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
)

