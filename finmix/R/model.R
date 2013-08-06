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

setClass("model",
         representation(dist        ="character",
		                r           = "integer",
		                K           = "integer",
		                weight      = "matrix",
		                par         = "list",
		                indicmod    = "character",
		                indicfix    = "logical",
		                T           = "matrix"),
         validity = function(object) {
             choices <- c("normal", "normult", "exponential", 
                          "student", "studmult", "poisson", 
                          "cond.poisson", "binomial")
             univ.choices <-  c("normal", "exponential", "student", 
                                "poisson", "cond.poisson", "binomial")	
			 multiv.choices <- c("normult", "studmult")
			 indicmod.choices <- c("multinomial")
			 if (!(object@dist %in% choices)) {
                 stop("Unknown distribution.")
             }
			 if (object@K <= 0) {					
                 stop("Number of components 'K' must be a positive integer.")
             }
			 if (object@indicmod != "multinomial") {
                 stop("Distribution of the indicator model 'indicmod' must be 
                      'multinomial'.")
             }
             if (object@K != ncol(object@weight)) {
                 stop("Dimension of slot 'weight' does not match number of 
                      components.")
             }
		     ## else: OK ##
			 TRUE
		 }
)

## Constructor for class 'model' ##
"model" <- function(dist = "normal", r = as.integer(1), 
                    K = as.integer(1), weight = matrix(), 
                    par = list(), indicmod = "multinomial", 
           			indicfix = FALSE, T = matrix(1)) {
    if(K > 1 && all(is.na(weight))) {
	    weight <- matrix(1/K, nrow = 1, ncol = K)
	}
	object <- new("model", dist = dist, r = as.integer(r), 
                 K = as.integer(K), weight = weight, par = par, 
				 indicmod = indicmod, indicfix = indicfix, T = T)
    return(object)
}

setMethod("plot", "model", function(x, y, ..., dis.grid = 1:10, persp.grid.x = seq(-10, 10, length = 40), 
					persp.grid.y = seq(-10, 10, length = 40), theta = 0, phi = 15, deparse.level = 1) {
					object <- x
					if(object@dist == "normal") {
						main <- paste("Normal Mixture K=", object@K, sep="")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:object@K) {
								fun.value <- fun.value + object@weight[i] * dnorm(y, 
									mean = object@par$mu[i], sd = object@par$sigma[i]) 
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(object@dist == "normult") {
						main.outer <- paste("Mutlivariate Normal Mixture K=", object@K, sep ="")
							
						x.axis <- persp.grid.x
						y.axis <- persp.grid.y
						z.axis <- x.axis %*% t(y.axis)
						fun <- function(y, Inputobject) {
							fun.value <- 0
							for(i in 1:Inputobject@K) {
								fun.value <- fun.value + Inputobject@weight[i] * 
									dmnorm(y, mean = Inputobject@par$mu[,i], 
										varcov = Inputobject@par$sigma[,,i])
							}
							fun.value
						}
				
						if(object@r == 2) {
							for(i in 1:length(x.axis)) {
								for(j in 1:length(y.axis)){
									input <- as.vector(c(x.axis[i], y.axis[j]))
									z.axis[i,j] <- fun(input, object)
								}
							}
							par(mfcol = c(2,2))
							persp(x.axis, y.axis, z.axis, theta = theta, phi = phi, ...)
							contour(x.axis, y.axis, z.axis, ...)
							image(x.axis, y.axis, z.axis, ...)
							mtext(main.outer, outer = TRUE)
						}
						else if(object@r > 2 && object@r < 6) {
							par(mfcol = c(ceiling(choose(object@r, 2)/3), 3))
							par(mar = c(2, 2, 2, 0.5), oma = c(2, 2, 2, 2))
							browser(expr = TRUE)	
							for(k in 1:(object@r - 1)) {
								for(l in (k + 1):object@r) {
									marmodel <- mixturemar(object, J = c(k, l))
									for(i in 1:length(x.axis)) {
						 				for(j in 1:length(y.axis)) {
											input <- as.vector(c(x.axis[i], 
													y.axis[j]))
											z.axis[i, j] <- fun(input, marmodel) 
										}
									}
									main <- paste("var ", k, " & var ", l, sep ="")
									contour(x.axis, y.axis, z.axis, main = main, ...)
								}
								remove(marmodel)
							} 
							mtext(main.outer, outer = TRUE)	
						}
						else 
							return("'plot()' is not implemented for model of dimension > 5. Use 'mixturemar()' to create margin models and then use 'plot' to show marginal distributions.")
					}
					else if(object@dist == "exponential") {
						main <- paste("Exponential Mixture K=", object@K, sep = "")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:object@K) {
								fun.value <- fun.value + object@weight[i] * dexp(y, rate = 
										object@par$lambda[i])
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(object@dist == "student") {
						main <- paste("Student-t Mixture K=", object@K, sep="")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:object@K) {
								fun.value <- fun.value + object@weight[i] * dstud(y, mu = object@par$mu[i], sigma = object@par$sigma[i], df = object@par$df[i])
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(object@dist == "studmult") {
						main.outer <- paste("Mutlivariate Student-t Mixture K=", object@K, sep ="")
						x.axis <- persp.grid.x
						y.axis <- persp.grid.y
						z.axis <- x.axis %*% t(y.axis)
						fun <- function(y, Inputobject) {
							fun.value <- 0
							for(i in 1:Inputobject@K) {
								fun.value <- fun.value + Inputobject@weight[i] * 
									dmt(y, mean = Inputobject@par$mu[,i], S =
										Inputobject@par$sigma[,,i], 
										df = Inputobject@par$df[i])
							}
							fun.value
						}
				
						if(object@r == 2) {
							for(i in 1:length(x.axis)) {
								for(j in 1:length(y.axis)){
									input <- as.vector(c(x.axis[i], y.axis[j]))
									z.axis[i,j] <- fun(input, object)
								}
							}

							par(mfcol = c(2,2), mar = c(2,2,0.5,0.5), oma = c(2, 2, 2, 1))
							persp(x.axis, y.axis, z.axis, theta = theta, phi = phi, ...)
							contour(x.axis, y.axis, z.axis, ...)
							image(x.axis, y.axis, z.axis, ...)
							mtext(main.outer, outer = TRUE)
						}
						else if(object@r > 2 && object@r < 6) {
							browser(expr = TRUE)
							par(mfcol = c(ceiling(choose(object@r, 2)/3), 3))
							par(mar = c(2, 2, 2, 0.5), oma = c(2, 2, 2, 2))
									
							for(k in 1:(object@r - 1)) {
								for(l in (k + 1):object@r) {
									marmodel <- mixturemar(object, J = c(k, l))
									for(i in 1:length(x.axis)) {
						 				for(j in 1:length(y.axis)) {
											input <- as.vector(c(x.axis[i], 
													y.axis[j]))
											z.axis[i, j] <- fun(input, marmodel) 
										}
									}
									xlab.name <- paste("var ", k, sep = "")
									ylab.name <- paste("var ", l, sep = "")
									main <- paste("var ", k, " & var ", l, sep ="")
									contour(x.axis, y.axis, z.axis, main = main, ...)
								}
								remove(marmodel)
							} 
							mtext(main.outer, outer = TRUE)	
						}
						else 
							return("'plot()' is not implemented for model of dimension > 5. Use 'mixturemar()' to create margin models and then use 'plot' to show marginal distributions.")
	
					}
					else if(object@dist == "poisson" || object@dist == "cond.poisson") {
						main <- paste("Poisson Mixture K=", object@K, sep="")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:object@K) {
								fun.value <- fun.value + object@weight[i] * dpois(y, 
										lambda = object@par$lambda[i])
							}
							fun.value
						}
						y <- dis.grid
						y.fun <- fun(y)
						names.grid <- as.character(dis.grid)
						barplot(y.fun, main = main, names.arg = names.grid, ...)
					}
					else if(object@dist == "binomial") {
						main <- paste("Binomial Mixture K=", object@K, sep="")
						if(length(object@T) != 1)
							return("[Error] Plotting a binomial distribution with differing repetitions is not possible.")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:object@K) {
								fun.value <- fun.value + object@weight[i] * dbinom(y, p =
									object@par$p[i], size = object@T[1])
							}
		
							return(fun.value)
						}
						y <- 1:object@T[1]
						y.fun <- fun(y)
						names.grid <- as.character(1:object@T[1])
						barplot(y.fun, main = main, names.arg = names.grid, ...)
					}
				}
)

## Marginal Mixture ##
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))
setMethod("mixturemar", "model", function(object, J) {
					if(object@dist == "normult") {	
						dist <- ifelse(length(J) == 1, "normal", "normult")
						r <- length(J)
						K <- object@K
						weight <- object@weight
						mu <- array(object@par$mu[J,1:K], dim = c(r, K))
						if(r == 1) {
							sigma <- array(object@par$sigma[J, J, 1:K], dim =c(r, K))
						}
						else {
							sigma <- array(object@par$sigma[J, J, 1:K], dim = c(r, r, K))
						}
						par <- list(mu = mu, sigma = sigma)
						indicmod <- "multinomial"
						indicfix <- TRUE
						T <- matrix()
						margin.model <- new("model", dist = dist, r = r, K = K, weight = weight, 
									par = par, indicmod = indicmod, 
									indicfix = indicfix, T = T)
						valiobject(margin.model)
						return(margin.model)
					}
					else if(object@dist == "studmult") {
						dist <- ifelse(length(J) == 1, "student", "studmult")
						r <- length(J)
						K <- object@K
						weight <- object@weight
						if(r == 1) {
							df <- object@par$df * matrix(object@par$sigma[J,J,1:K], nrow = 1,
											 ncol = K)
						}
						else {
							df <- object@par$df
						}
						mu <- array(object@par$mu[J,1:K], dim = c(r, K))
						if(r > 1) {
							sigma <- array(object@par$sigma[J, J, 1:K], dim =c(r, r, K))
						}
						if(r == 1) {
							par <- list(df = df, ncp = mu)
						}
						else  {
							par <- list(df = df, mu = mu, sigma = sigma)
						}
						indicmod <- "multinomial"
						indicfix <- TRUE
						T <- matrix()
						margin.model <- new("model", dist = dist, r = r, K = K, weight = weight, 
									par = par, indicmod = indicmod,
									indicfix = indicfix, T = T)
						valiobject(margin.model)
						return(margin.model)

					}
					else {
						return("[Error]: Marginal distribution can only be obtained from multivariate distribution.")
					}
				}
)
	
## Show ##
setMethod("show", "model", 
          function(object) {
              cat("Object 'model'\n")
              cat("     class       :", class(object), "\n")
              cat("     dist        :", object@dist, "\n")
              cat("     r           :", object@r, "\n")
              cat("     K           :", object@K, "\n")
              cat("     weight      :", 
                  paste(dim(object@weight), collapse = "x"), "\n")
              if (!all(is.na(object@par))) {
                  cat("     par         : List of",
                      length(object@par), ",\n")
              }
              cat("     indicmod    :", object@indicmod, "\n")
              cat("     indicfix    :", object@indicfix, "\n")
              if (object@dist == 'binomial') {
                  if(length(object@T) > 1) {
                      cat("     T       :",
                          paste(dim(object@T), collapse = "x"), "\n")
                  } else {
                      cat("     T       :", object@T, "\n")
                  }
              }				
          }
)
## Getters ##
setGeneric("getDist", function(object) standardGeneric("getDist"))
setMethod("getDist", "model", function(object) {
					return(object@dist)
				}
)
setGeneric("getR", function(object) standardGeneric("getR"))
setMethod("getR", "model", function(object) {
					return(object@r)
				}
)
setGeneric("getK", function(object) standardGeneric("getK"))
setMethod("getK", "model", function(object) {
					return(object@K)
				}
)
setGeneric("getWeight", function(object) standardGeneric("getWeight"))
setMethod("getWeight", "model", function(object) {
						return(object@weight)
					}
)
setGeneric("getPar", function(object) standardGeneric("getPar"))
setMethod("getPar", "model", function(object) {
					return(object@par)
				}
)
setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))
setMethod("getIndicmod", "model", function(object) {
						return(object@indicmod)					
					}
)
setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))
setMethod("getIndicfix", "model", function(object) {
						return(object@indicfix)
					}
)
setGeneric("getT", function(object) standardGeneric("getT"))
setMethod("getT", "model", function(object) {
					return(object@T)
				}
)
## Setters ##
setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))
setReplaceMethod("setDist", "model", function(object, value) {
						object@dist <- value
						validObject(object)
						return(object)
					}
)
setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))
setReplaceMethod("setR", "model", function(object, value){
						object@r <- as.integer(value)
						validObject(object)
						return(object)
					}
)
setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))
setReplaceMethod("setK", "model", function(object, value) {
						object@K <- as.integer(value)
						object@weight <- matrix(1/value, 
							nrow = 1, ncol = value)
						validObject(object)
						return(object)
					}
)
setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))
setReplaceMethod("setWeight", "model", function(object, value) {
						object@weight <- value
						validObject(object)
						return(object)
					}
)
setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))
setReplaceMethod("setPar", "model", function(object, value) {
						object@par <- value
						validObject(object)
						return(object)
					}
)
setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))
setReplaceMethod("setIndicmod", "model", function(object, value) {
					object@indicmod <- value
					validObject(object)
					return(object)
				}
)
setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))
setReplaceMethod("setIndicfix", "model", function(object, value) {
							object@indicfix <- value
							validObject(object)
							return(object)
						}
)
setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))
setReplaceMethod("setT", "model", function(object, value) {
						object@T <- value
						validObject(object)
						return(object)
					}
)
