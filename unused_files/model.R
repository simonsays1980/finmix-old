setClass("model",
	representation(
		dist="character",
		r = "numeric",
		K = "numeric",
		weight = "matrix",
		par = "list",
		indicmod = "numeric"
	),
	validity = function(object) {
				choices <- c("normal", "normult", "exponential", "student", "studmult", "poisson", "binomial")
				is.dist <- object@dist %in% choices
				if(!is.dist) 
					return("Unknown distribution.")
				if(object@r == 1 && (object@dist %in% c("normult", "studmult")))
					return("Univariate dimension but multivariate distribution.")
				univ.choices <-  c("normal", "exponential", "student", "poisson", "binomial")	
				if(object@r > 1 && (object@dist %in% univ.choices)) 
					return("Multivariate dimension but univariate distribution.")
				# else: OK ##

				TRUE
			}
)

## Constructor for class 'model' ##
"model" <- function(dist = "normal", r = 1, K = 1, weight = matrix(), par = list()) {
		model <- new("model", dist = dist, r = r, K = K, weight = weight, par = par)
		return(model)
}

## Initialize function ##
setMethod("initialize", "model", function(.Object, dist = "normal", r = 1, K = 1, weight = matrix(), par = list(NA)) {
						## set the parameter and weight matrices for the choices ##
						weight.def <- matrix(1/K, nrow = 1, ncol = K)
						if(dist == "normal") {
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								.Object@par <- list(mu = array(0, dim = c(1, K)),
									sigma = array(1, dim = c(1, K)))
							}
							else {
								.Object@par <- par
								names(.Object@par) <- c("mu", "sigma")
							}
						}
						else if(dist == "normult") {
							r <- ifelse(r == 1, 2, r)
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								mu.def <- array(0, dim = c(r, K))
								sigma.def.cov <- array(diag(1, ncol = r, nrow = r))
								.Object@par <- list(mu = mu.def,
									 sigma = array(sigma.def.cov, dim = c(r, r, K)))
							}
							else {
								.Object@par <- par
								names(.Object@par) <- c("mu", "sigma")
							}
						}
						else if(dist == "exponential") {
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								.Object@par <- list(lambda = array(1, dim = c(1, K)))
							}
							else {
								.Object@par <- par
								names(.Object@par) <- c("lambda")
							}
						}
						else if(dist == "student") {
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								.Object@par <- list(mu = array(0, dim = c(1, K)), 
											sigma = array(1, dim = c(1, K)),
											df = array(100, dim = c(1, K)))
							}
							else {
								.Object@par <- par 
								names(.Object@par) <- c("mu", "sigma", "df")
							}
						}
						else if(dist == "studmult") {
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								df.def <- array(Inf, dim = c(1, K))
								mu.def <- array(0, dim=c(r, K))
								sigma.def <- array(diag(1, nrow = r, ncol = r), dim =c(r,r,K))
								.Object@par <- list(df = df.def, mu = mu.def, 
											sigma = sigma.def)
							}
							else {
								.Object@par <- par
								names(.Object@par) <- c("df", "mu", "sigma")	
							}
						}
						else if(dist == "poisson") { 
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								.Object@par <- list(lambda = array(1, dim = c(1, K)))
							}
							else {
								.Object@par <- par
								names(.Object@par) <- c("lambda")
							}
						}
						else if(dist == "binomial") {
							if(all(is.na(weight))) {
								.Object@weight <- weight.def
							}	
							else {
								.Object@weight <- weight
							}
							if(all(is.na(par))) {
								.Object@par <- list(p = array(0.5, dim = c(1, K)),
										size = array(10, dim = c(1, K)))
							}
							else {
								.Object@par <- par
								names(.Object@par) <- c("n", "p")
							}
						}
						.Object@dist <- dist
						.Object@K <- K
						.Object@r <- r
						callNextMethod(.Object)
						return(.Object)
					}
)

setMethod("plot", "model", function(x, ..., dis.grid = 1:10, persp.grid = seq(-2, 2, length = 40), deparse.level = 1) {
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
						x.axis <- persp.grid
						y.axis <- persp.grid
						z.axis <- x.axis %*% t(y.axis)
						fun <- function(y, Input.Object) {
							fun.value <- 0
							for(i in 1:Input.Object@K) {
								fun.value <- fun.value + Input.Object@weight[i] * 
									dmnorm(y, mean = Input.Object@par$mu[,i], varcov =
									Input.Object@par$sigma[,,i])
							}
							fun.value
						}
						for(i in 1:length(x.axis)) {
							for(j in 1:length(y.axis)){
								input <- as.vector(c(x.axis[i], y.axis[j]))
								z.axis[i,j] <- fun(input, .Object)
							}
						}

						if(.Object@r == 2) {
							par(mfcol = c(2,2))
							persp(x.axis, y.axis, z.axis, ...)
							contour(x.axis, y.axis, z.axis, ...)
							image(x.axis, y.axis, z.axis, ...)
							mtext(main.outer, outer = TRUE)
						}
						else if(.Object@r > 2 && .Object@r < 6) {
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
								fun.value <- fun.value + .Object@weight[i] * dt(y, df = 
									.Object@par$df[i], ncp = .Object@par$ncp[i])
							}
							fun.value
						}
						curve(fun, main = main, ...)
					}
					else if(.Object@dist == "studmult") {
						main.outer <- paste("Mutlivariate Student-t Mixture K=", .Object@K, sep ="")
						x.axis <- persp.grid
						y.axis <- persp.grid
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
						for(i in 1:length(x.axis)) {
							for(j in 1:length(y.axis)){
								input <- as.vector(c(x.axis[i], y.axis[j]))
								z.axis[i,j] <- fun(input, .Object)
							}
						}

						if(.Object@r == 2) {
							par(mfcol = c(2,2), mar = c(2,2,0.5,0.5), oma = c(2, 2, 2, 1))
							persp(x.axis, y.axis, z.axis, ...)
							contour(x.axis, y.axis, z.axis, ...)
							image(x.axis, y.axis, z.axis, ...)
							mtext(main.outer, outer = TRUE)
						}
						else if(.Object@r > 2 && .Object@r < 6) {
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
					else if(.Object@dist == "poisson") {
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
						fun <- function(y) {
							fun.value <- 0
							for(i in 1:.Object@K) {
								fun.value <- fun.value + .Object@weight[i] * dbinom(y, p =
									.Object@par$p[i], size = .Object@par$size[i])
							}
							fun.value
						}
						y <- dis.grid
						y.fun <- fun(y)
						names.grid <- as.character(dis.grid)
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
						margin.model <- new("model", dist = dist, r = r, K = K, weight = weight, 
									par = par)
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
						margin.model <- new("model", dist = dist, r = r, K = K, weight = weight, 
									par = par)
						validObject(margin.model)
						return(margin.model)

					}
					else {
						return("[Error]: Marginal distribution can only be obtained from multivariate distribution.")
					}
				}
)
	
## Show ##
setMethod("show", "model", function(.Object) {
					cat("Object 'model'\n")
					cat("	Class	: ", class(.Object), "\n")
					cat("	Dist 	: ", .Object@dist, "\n")
					cat("	r	: ", .Object@r, "\n")
					cat("	K	: ", .Object@K, "\n")
					cat("	weight	:  [", .Object@weight, "]\n")
					cat("	Par	:  List of ", length(names(getPar(.Object))), "\n")
					if(.Object@dist == "normal") {
						cat("	mu	: [", .Object@par$mu, "]\n")
						cat("	sigma	: [", .Object@par$sigma, "]\n")
					}
					else if(.Object@dist == "normult") {
						cat("	mu	: [", .Object@par$mu[1,], "]\n")
						for(i in 2:.Object@r) {
							cat("		  [", .Object@par$mu[i,], "]\n")
						}
						cat("	sigma	: K = 1\n")
						for(i in 1:.Object@K) {
							if(i != 1) {
								cat("		  K = ", i, "\n")
							}
							for(j in 1:.Object@r){
								cat("		  [", .Object@par$sigma[j,,i], "]\n")
							}
						}
					}
					else if(.Object@dist == "exponential") {
						cat("	lambda	: [", .Object@par$lambda, "]\n")
					}
					else if(.Object@dist == "student") {
						cat("	df	: [", .Object@par$df, "]\n")
						cat("	ncp	: [", .Object@par$ncp, "]\n")
					}
					else if(.Object@dist == "studmult") {
						cat("	df	: [", .Object@par$df, "]\n")
						cat("	mu	: [", .Object@par$mu[1, ], "]\n")
						for(i in 2:.Object@r) {
							cat("		  [", .Object@par$mu[i, ], "]\n")
						}
						cat("	sigma 	: K = 1\n")
						for(i in 1:.Object@K) {
							if(i != 1) {
								cat("		  K = ",i,"\n", sep = "")
							}
							for(j in 1:.Object@r) {
								cat("		   [", .Object@par$sigma[j,,i], "]\n")
							}
						}
					}
					else if(.Object@dist == "poisson") {
						cat("	lambda	: [", .Object@par$lambda, "]\n")
					}
					else if(.Object@dist == "binomial") {
						cat("	p	: [", .Object@par$p, "]\n")
						cat("	size	: [", .Object@par$size, "]\n")
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
						return(.Objec@weight)
					}
)
setGeneric("getPar", function(.Object) standardGeneric("getPar"))
setMethod("getPar", "model", function(.Object) {
					return(.Object@par)
				}
)

## Setters ##
setGeneric("setDist<-", function(.Object, value) standardGeneric("setDist<-"))
setReplaceMethod("setDist", "model", function(.Object, value) {
						.Object@dist <- value[1]
						.Object@r <- value[2]
						validObject(.Object)
						return(.Object)
					}
)
## ( OFF ) This setter is inactivated due to a more secure handling by users ##
## setGeneric("setR<-", function(.Object, value) standardGeneric("setR<-"))
## setReplaceMethod("setR", "model", function(.Object, value){
##						.Object@r <- value
##						validObject(.Object)
##						return(.Object)
##					}
##)
setGeneric("setK<-", function(.Object, value) standardGeneric("setK<-"))
setReplaceMethod("setK", "model", function(.Object, value) {
						.Object@K <- value
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
