setClass("model",
	representation(
		dist="character",
		r = "numeric",
		K = "numeric",
		weight = "matrix",
		par = "list",
		indicmod = "character",
		indicfix = "logical",
		T = "matrix",
		.cache = "environment"
	),
	validity = function(object) {
				choices <- c("normal", "normult", "exponential", "student", 
					"studmult", "poisson", "cond.poisson", "binomial")
				univ.choices <-  c("normal", "exponential", "student", 
					"poisson", "cond.poisson", "binomial")	
				multiv.choices <- c("normult", "studmult")
				indicmod.choices <- c("multinomial")
				if(!(object@dist %in% choices)) 
					return("[Error] Unknown distribution.")
				if(object@K <= 0) 
					return("[Error] Number of components 'K' must be a positive integer.")
				if(object@indicmod != "multinomial")
					return("[Error] Distribution of the indicator model 'indicmod' must be Multinomial.")
				if(object@K != ncol(object@weight))
					return("[Error] Dimension of weight does not match number of components.")
				# else: OK ##

				TRUE
			}
)

## Initializer ##
setMethod("initialize", "model", function(.Object, ..., .cache = new.env()){
					callNextMethod(.Object, .cache = .cache, ...)				
				}
)

## Constructor for class 'model' ##
"model" <- function(dist. = "normal", r. = 1, K. = 1, weight. = matrix(), par. = list(), indicmod. = "multinomial", 
			indicfix. = FALSE, T. = matrix(1)) {

				if(K. > 1 && all(is.na(weight.))) {
					weight. <- matrix(1/K., nrow = 1, ncol = K.)
				}
				model <- new("model", dist = dist., r = r., K = K., weight = weight., par = par., 
					indicmod = indicmod.,indicfix = indicfix., T = T.)
		return(model)
}

setMethod("plot", "model", function(x, y, ..., dis.grid = 1:10, persp.grid.x = seq(-10, 10, length = 40), 
					persp.grid.y = seq(-10, 10, length = 40), theta = 0, phi = 15, deparse.level = 1) {
					.Object <- x
					if(.Object@dist == "normal") {
						main <- paste("Normal Mixture K=", .Object@K, sep="")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:.Object@K) {
								fun.value <- fun.value + .Object@weight[i] * dnorm(y, 
									mean = .Object@par$mu[i], sd = .Object@par$sigma[i]) 
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(.Object@dist == "normult") {
						main.outer <- paste("Mutlivariate Normal Mixture K=", .Object@K, sep ="")
							
						x.axis <- persp.grid.x
						y.axis <- persp.grid.y
						z.axis <- x.axis %*% t(y.axis)
						fun <- function(y, Input.Object) {
							fun.value <- 0
							for(i in 1:Input.Object@K) {
								fun.value <- fun.value + Input.Object@weight[i] * 
									dmnorm(y, mean = Input.Object@par$mu[,i], 
										varcov = Input.Object@par$sigma[,,i])
							}
							fun.value
						}
				
						if(.Object@r == 2) {
							for(i in 1:length(x.axis)) {
								for(j in 1:length(y.axis)){
									input <- as.vector(c(x.axis[i], y.axis[j]))
									z.axis[i,j] <- fun(input, .Object)
								}
							}
							par(mfcol = c(2,2))
							persp(x.axis, y.axis, z.axis, theta = theta, phi = phi, ...)
							contour(x.axis, y.axis, z.axis, ...)
							image(x.axis, y.axis, z.axis, ...)
							mtext(main.outer, outer = TRUE)
						}
						else if(.Object@r > 2 && .Object@r < 6) {
							par(mfcol = c(ceiling(choose(.Object@r, 2)/3), 3))
							par(mar = c(2, 2, 2, 0.5), oma = c(2, 2, 2, 2))
							browser(expr = TRUE)	
							for(k in 1:(.Object@r - 1)) {
								for(l in (k + 1):.Object@r) {
									marmodel <- mixturemar(.Object, J = c(k, l))
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
					else if(.Object@dist == "exponential") {
						main <- paste("Exponential Mixture K=", .Object@K, sep = "")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:.Object@K) {
								fun.value <- fun.value + .Object@weight[i] * dexp(y, rate = 
										.Object@par$lambda[i])
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(.Object@dist == "student") {
						main <- paste("Student-t Mixture K=", .Object@K, sep="")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:.Object@K) {
								fun.value <- fun.value + .Object@weight[i] * dstud(y, mu = .Object@par$mu[i], sigma = .Object@par$sigma[i], df = .Object@par$df[i])
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(.Object@dist == "studmult") {
						main.outer <- paste("Mutlivariate Student-t Mixture K=", .Object@K, sep ="")
						x.axis <- persp.grid.x
						y.axis <- persp.grid.y
						z.axis <- x.axis %*% t(y.axis)
						fun <- function(y, Input.Object) {
							fun.value <- 0
							for(i in 1:Input.Object@K) {
								fun.value <- fun.value + Input.Object@weight[i] * 
									dmt(y, mean = Input.Object@par$mu[,i], S =
										Input.Object@par$sigma[,,i], 
										df = Input.Object@par$df[i])
							}
							fun.value
						}
				
						if(.Object@r == 2) {
							for(i in 1:length(x.axis)) {
								for(j in 1:length(y.axis)){
									input <- as.vector(c(x.axis[i], y.axis[j]))
									z.axis[i,j] <- fun(input, .Object)
								}
							}

							par(mfcol = c(2,2), mar = c(2,2,0.5,0.5), oma = c(2, 2, 2, 1))
							persp(x.axis, y.axis, z.axis, theta = theta, phi = phi, ...)
							contour(x.axis, y.axis, z.axis, ...)
							image(x.axis, y.axis, z.axis, ...)
							mtext(main.outer, outer = TRUE)
						}
						else if(.Object@r > 2 && .Object@r < 6) {
							browser(expr = TRUE)
							par(mfcol = c(ceiling(choose(.Object@r, 2)/3), 3))
							par(mar = c(2, 2, 2, 0.5), oma = c(2, 2, 2, 2))
									
							for(k in 1:(.Object@r - 1)) {
								for(l in (k + 1):.Object@r) {
									marmodel <- mixturemar(.Object, J = c(k, l))
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
					else if(.Object@dist == "poisson" || .Object@dist == "cond.poisson") {
						main <- paste("Poisson Mixture K=", .Object@K, sep="")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:.Object@K) {
								fun.value <- fun.value + .Object@weight[i] * dpois(y, 
										lambda = .Object@par$lambda[i])
							}
							fun.value
						}
						y <- dis.grid
						y.fun <- fun(y)
						names.grid <- as.character(dis.grid)
						barplot(y.fun, main = main, names.arg = names.grid, ...)
					}
					else if(.Object@dist == "binomial") {
						main <- paste("Binomial Mixture K=", .Object@K, sep="")
						if(length(.Object@T) != 1)
							return("[Error] Plotting a binomial distribution with differing repetitions is not possible.")
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:.Object@K) {
								fun.value <- fun.value + .Object@weight[i] * dbinom(y, p =
									.Object@par$p[i], size = .Object@T[1])
							}
		
							return(fun.value)
						}
						y <- 1:.Object@T[1]
						y.fun <- fun(y)
						names.grid <- as.character(1:.Object@T[1])
						barplot(y.fun, main = main, names.arg = names.grid, ...)
					}
				}
)

## Marginal Mixture ##
setGeneric("mixturemar", function(.Object, J) standardGeneric("mixturemar"))
setMethod("mixturemar", "model", function(.Object, J) {
					if(.Object@dist == "normult") {	
						dist <- ifelse(length(J) == 1, "normal", "normult")
						r <- length(J)
						K <- .Object@K
						weight <- .Object@weight
						mu <- array(.Object@par$mu[J,1:K], dim = c(r, K))
						if(r == 1) {
							sigma <- array(.Object@par$sigma[J, J, 1:K], dim =c(r, K))
						}
						else {
							sigma <- array(.Object@par$sigma[J, J, 1:K], dim = c(r, r, K))
						}
						par <- list(mu = mu, sigma = sigma)
						indicmod <- "multinomial"
						indicfix <- TRUE
						T <- matrix()
						margin.model <- new("model", dist = dist, r = r, K = K, weight = weight, 
									par = par, indicmod = indicmod, 
									indicfix = indicfix, T = T)
						validObject(margin.model)
						return(margin.model)
					}
					else if(.Object@dist == "studmult") {
						dist <- ifelse(length(J) == 1, "student", "studmult")
						r <- length(J)
						K <- .Object@K
						weight <- .Object@weight
						if(r == 1) {
							df <- .Object@par$df * matrix(.Object@par$sigma[J,J,1:K], nrow = 1,
											 ncol = K)
						}
						else {
							df <- .Object@par$df
						}
						mu <- array(.Object@par$mu[J,1:K], dim = c(r, K))
						if(r > 1) {
							sigma <- array(.Object@par$sigma[J, J, 1:K], dim =c(r, r, K))
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
						validObject(margin.model)
						return(margin.model)

					}
					else {
						return("[Error]: Marginal distribution can only be obtained from multivariate distribution.")
					}
				}
)
	
## Show ##
setMethod("show", "model", function(object) {
					cat("Object 'model'\n")
					cat("	type	    :", class(object), "\n")
					cat("	dist 	    :", object@dist, "\n")
					cat("	r	        :", object@r, "\n")
					cat("	K	        :", object@K, "\n")
					cat("	weight	    :", paste(dim(object@weight), collapse = "x"), "\n")
					if(!all(is.na(object@par))) {
						cat("	par	        : List of ", length(names(getPar(object))), "\n")
					}	
					cat("	indicmod    :", object@indicmod,"\n")
					cat("	indicfix    :", object@indicfix, "\n")
					if(object@dist == "binomial") {
						if(length(T) > 1) {
							cat("	T	        :", paste(dim(object@T), collapse = "x"), "\n")
						}
						else { ## T is the same for all 
							cat("	T	        :", object@T, "\n")
						}
					}
				}
)
## Getters ##
setGeneric("getDist", function(.Object) standardGeneric("getDist"))
setMethod("getDist", "model", function(.Object) {
					return(.Object@dist)
				}
)
setGeneric("getR", function(.Object) standardGeneric("getR"))
setMethod("getR", "model", function(.Object) {
					return(.Object@r)
				}
)
setGeneric("getK", function(.Object) standardGeneric("getK"))
setMethod("getK", "model", function(.Object) {
					return(.Object@K)
				}
)
setGeneric("getWeight", function(.Object) standardGeneric("getWeight"))
setMethod("getWeight", "model", function(.Object) {
						return(.Object@weight)
					}
)
setGeneric("getPar", function(.Object) standardGeneric("getPar"))
setMethod("getPar", "model", function(.Object) {
					return(.Object@par)
				}
)
setGeneric("getIndicMod", function(.Object) standardGeneric("getIndicMod"))
setMethod("getIndicMod", "model", function(.Object) {
						return(.Object@indicmod)					
					}
)
setGeneric("getIndicFix", function(.Object) standardGeneric("getIndicFix"))
setMethod("getIndicFix", "model", function(.Object) {
						return(.Object@indicfix)
					}
)
setGeneric("getT", function(.Object) standardGeneric("getT"))
setMethod("getT", "model", function(.Object) {
					return(.Object@T)
				}
)
## Setters ##
setGeneric("setDist<-", function(.Object, value) standardGeneric("setDist<-"))
setReplaceMethod("setDist", "model", function(.Object, value) {
						.Object@dist <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setR<-", function(.Object, value) standardGeneric("setR<-"))
setReplaceMethod("setR", "model", function(.Object, value){
						.Object@r <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setK<-", function(.Object, value) standardGeneric("setK<-"))
setReplaceMethod("setK", "model", function(.Object, value) {
						.Object@K <- value
						.Object@weight <- matrix(1/value, 
							nrow = 1, ncol = value)
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setWeight<-", function(.Object, value) standardGeneric("setWeight<-"))
setReplaceMethod("setWeight", "model", function(.Object, value) {
						.Object@weight <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setPar<-", function(.Object, value) standardGeneric("setPar<-"))
setReplaceMethod("setPar", "model", function(.Object, value) {
						.Object@par <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setIndicMod<-", function(.Object, value) standardGeneric("setIndicMod<-"))
setReplaceMethod("setIndicMod", "model", function(.Object, value) {
					.Object@indicmod <- value
					validObject(.Object)
					return(.Object)
				}
)
setGeneric("setIndicFix<-", function(.Object, value) standardGeneric("setIndicFix<-"))
setReplaceMethod("setIndicFix", "model", function(.Object, value) {
							.Object@indicfix <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setT<-", function(.Object, value) standardGeneric("setT<-"))
setReplaceMethod("setT", "model", function(.Object, value) {
						.Object@T <- value
						validObject(.Object)
						return(.Object)
					}
)
