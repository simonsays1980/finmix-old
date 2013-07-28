setClass("prior",
	representation(
	weight 	= "matrix",
	par 	= "list",
	type 	= "character",
	hier 	= "logical"
	),
	validity = function(object) {
			type.choices <- c("condconjugate", "independent")
			if(!(object@type %in% type.choices)) 
				stop("Unknown prior 'type'. 'type' must be 'independent' 
                     or 'condconjugate'.")
			## else: OK
			TRUE		
	}
)

"prior" <- function(weight = matrix(), par = list(), type = "independent", 
                    hier = TRUE) {
			prior <- new("prior", weight = weight, par = par, type = type, 
                         hier = hier)
			return(prior)
}
"priordefine" <- function(data = data(), model = model(), coef.mat = NULL, 
                          varargin = NULL) {
	model.choices <- c("normal", "normult", "exponential", "student",
		"studmult", "poisson", "cond.poisson", "binomial")
	if (!(model@dist %in% model.choices))
		stop("The field 'dist' is obligatory'. It must be one of 'normal', 
               'normult', 'exponential', 'student', 'studmult', 'poisson', 
               'cond.poisson'  or 'binomial'.")
	if (model@dist == "cond.poisson") {
		if (is.null(coef.mat)) {
			stop("For a conditional Poisson mixture a coefficient matrix 
                 'coef.mat' has to be provided.") 
        } else if (!is.null(coef.mat)) {
            if (!is.matrix(coef.mat) && !is.array(coef.mat)) {
    			stop("Argument 'coef.mat' must be of type 'matrix' or 'array'.")
            } else if (nrow(coef.mat) != ncol(coef.mat)) {
    			stop("Argument 'coef.mat' must be a quadratic 'matrix' or 'array'.")
            } else if (nrow(coef.mat) != model@K || ncol(coef.mat) != model@K) {
    			stop("Dimension of argument 'coef.mat' must correspond to number 
                     of components 'K' in 'model'.\n")
            } else if (!(all(diag(coef.mat) == 1))) {
    			stop("Coefficients on the diagonal of 'coef.mat' must be equal 
                     to one.\n")
            }
        }
	}
	dist <- getDist(model)
	K <- getK(model)
	r <- getR(data)
	hier <- TRUE
	type <- "independent"
		if(!data@bycolumn) {
			datam <- t(data@y)
		}
		else {
			datam <- data@y
		}
		
		## poisson mixture ##
		if(dist == "poisson") {
			if(is.null(varargin)) {
				hier <- TRUE ## default prior is hierarchical
				type <- "condconjugate"
			}
			else {
				## check if object is valid ##
				if(!(class(varargin) == "prior")) { 
					stop("'varargin' must be of class 'prior'.")
                }
				validObject(varargin)
				hier <- varargin@hier
				type <- "condconjugate"
			}
			
			## default prior based on matching moments ##

			## choose level of overdispersion, depending on the ratio overdispersion/mean^2 ##
			## no idea data-based choice 
			mean <- mean(datam, na.rm = TRUE) 
			over <- var(datam, na.rm=TRUE) - mean
			if(over > 0) {
				a0 <- mean^2/over
			}
			else {
				a0 <- 10
			}
			
			if(hier) {
				g0 <- 0.5
				G0 <- mean * g0/a0
				be <- g0/G0
				par <- list(a = array(a0, dim = c(1, K)),
					b = array(be, dim = c(1, K)),
					g = g0, G = G0)
			}
			else {
				be = a0/mean
				par <- list(a = array(a0, dim = c(1, K)),
					b = array(be, dim = c(1, K)))
			}
		}
		
		## conditional poisson mixture ##
		else if (dist == "cond.poisson") {
			if (is.null(varargin)) {
				hier <- TRUE ## default prior is hierarchical
				type <- "condconjugate"
			} else {

				## check if object is valid ##
				if (!(class(varargin) == "prior")) {
					stop("'varargin' must be of class 'prior'.")
				}
				validObject(varargin)
				hier <- varargin@hier
				type <- varargin@hier 
			}
			if (type == "trunc.norm") {
				s0 <- sd(datam)
				par <- list(s = array(s0, dim = c(1,K)))
			} else if (type == "condconjugate") {
				mean <- mean(datam, na.rm = TRUE)
				over <- var(datam, na.rm = TRUE) - mean
				if (over > 0) {
					a0 <- mean^2/over
				} else {
					a0 <- 10
				}
				if(hier) {
					g0 <- 0.5
					G0 <- mean * g0/a0
					be <- g0/G0
					par <- list(a = array(a0, dim = c(1, K)),
						b = array(be, dim = c(1, K)),
						g = g0, G = G0)
					par$a <- par$a / 2^(rowSums(coef.mat) - 1)
					par$coef <- coef.mat
				} else {
					be <- a0/mean
					par <- list(a = array(a0, dim = c(1, K)),
						b = array(be, dim = c(1, K)))
					par$a <- par$a / 2^(rowSums(coef.mat) - 1)
					par$coef <- coef.mat
				}
			} 
		}
		## binomial mixture ##
		else if(dist == "binomial") {
			type <- "condconjugate"
			## uniform prior ##
			a0 <- 1
			b0 <- 1
			par <- list(a = array(a0, dim = c(1, K)),
				b = array(b0, dim = c(1, K)))
		}	
		## exponential mixture ##
		else if(dist == "exponential") {
			## prior following Wagner (2007) ##
			hier <- FALSE
			type <- "condconjugate"
			a0 <- 0.1
			be <- mean(datam) * a0
			par <- list(a = array(a0, dim = c(1, K)), 
				b = array(be, dim = c(1, K))) 
		}
		## normal or student-t mixtures ##
		else {
			## check if varargin is non-empty and prior object ##
			## set hierarchical or non-hierarchical prior ##
			if(is.null(varargin)) { 
				## default prior: independent hierarchical prior ##
				hier <- TRUE
				type <- "independent"
			}
			else {
				if(!(class(varargin) == "prior")) {
					return("[Error] 'varargin' must be of class 'prior'.")
				}
				else {
					validObject(varargin)
					hier <- ifelse(length(varargin@hier) > 0, varargin@hier, TRUE)
					type <- ifelse(length(varargin@type) > 0, varargin@type, TRUE)
				}
			}
 			conjugate.prior <- type == "condconjugate"
			bensmail <- FALSE      
			rich.green <- FALSE
			if(conjugate.prior || !hier) {
				bensmail <- TRUE ## Prior following Bensmail et al. 
			}
			else {
				rich.green <- TRUE ## Prior following Richardson and Green for r = 1
						   ##                 Stephens for r = 2 only 
			}
			
			if(rich.green) {
				## row vectors: dimension 1 x r
				max <- apply(datam, 2, max, na.rm = TRUE)
				min <- apply(datam, 2, min, na.rm = TRUE)
				mean <- (max + min) * 0.5
				cov <- diag((max - min)^2)
			}
			else {
				## row vectors: dimension 1 x r
				mean <- apply(datam, 2, mean, na.rm = TRUE)
				cov <- ifelse(r > 1, cov(datam), var(datam))
			}		
			b0 <- mean
			
			if(conjugate.prior) {
				B0sc <- 1 ## info contained in a standard conjugate (sc) prior (equal to N0)
			}
			else {
				B0inv <- solve(cov) ## info contained in a non-conjugate prior,
						    ## i.e. either by Richardson Green or by Benmail et al. 
			}
			if(!conjugate.prior) {
				if(r > 1) {
					par.mu <- list(b = array(t(b0), dim = c(r, K)),
							Binv = array(B0inv, dim = c(r, r, K)))
				}
				else { ## r = 1
					par.mu <- list(b = array(b0, dim = c(1,K)),
							Binv = array(B0inv, dim = c(1, K)))	
				}
			}
			else { ## conditionally conjugate prior
				if(r > 1) {
					par.mu <- list(b = array(t(b0), dim = c(r, K)), 
						N0 = array(B0sc, dim = c(r, r, K)))
				}
				else { ## r = 1
					par.mu <- list(b = array(b0, dim = c(1, K)),
							N0 = array(B0sc, dim = c(1, K)))
				}
			}
	
			## prior sigma ##
			## r = 1:	Inverse Gamma with c0, C0
			## r > 1:	Wishart with c0, C0
			## any r:	Q in {Inverse Gamma, Inverse Wishart} with prQnu (prior Q nu) and prQS
			##		We use the Gamma and Wishart and sample the inverse Variance.
			## where:	prQnu	degrees of freedom for Wishart and shape for Gamma
			## 	:	prQS	shape for Q
			## 	
			## Select Q0 the prior mean of Q. 
			## Determine prQS from prQS = Q0 * (prQnu - (r + 1)/2). This matches Q0 to the mean 
			## of the Inverse Gamma or the Inverse Wishart distribution and to the mode of Q {-1}
			## i.e. the Gamma and Wishart distribution respectively. 
			## Further, variance shrinkage towards the ratio prQS/dfQpr, where dfQpr bounds the 
			## ratio of the variances.
		
			dfQpr <- 2.5             ## this bounds the ratio of variances to 10 for r = 1
			prQnu <- dfQpr + (r-1)/2
			
			if(K == 1) {
				phi <- 1 ## c0 heterogeneity 
			}
			else {
				## Tuning of the prior for sigma is done by explained heterogeneity
				## See p. 192, chapter 6.3.2 Fruewirth-Schnatter (2006)
				## Rhet 
				## -> 1: 	means very different in relation to variances
				## -> 0: 	means rather similar in relation to variances 
				## 0 < Rhet < 1 (do not choose 0 nor 1)
				## SMALL VALUE: 	leads to very informative prior for mu_k
				##			close to b0. Should be chosen only in 
				##			combination with a hierarchical prior
				##			on b0.
				## LARGE VALUE:		leads to a very informative prior for 
				##			sigma_k close to prQS/prQnu. Should only
				##			be chosen in combination with hierarchical
				##			prior on prQS.
				Rhet <- 0.5 ## Rhet = 2/3
				phi <- (1 - Rhet) 
			}
			prQS <- cov * phi * (prQnu - (r + 1)/2) 
			
			if(r > 1) {
				detprQS <- log(det(matrix(prQS)))
			}
				
			if(hier) {
				if(rich.green) {
					if(r == 1) {
						g0 <- 0.2  ## Richardson and Green. Sampling from Gamma allows
							   ## arbitrary g0: 
							   ## WARNING: seems to cause problems in bayesf
						prQnu <- 2 ## Note that prQnu standard is changed here
					}
					else if(r == 2){
						g0 <- 0.3  ## Stephens	
							   ## WARNING:  seems to cause problems in bayesf
						prQnu <- 3 ## prQnu is changed also in relation from standard
					}
					else { ## r > 2
						g0 <- 0.5 + (r - 1)/2
					}
					g0 <- 0.5 + (r - 1)/2
					G0 <- 100 * g0/prQnu * solve(cov) ## Stephens
					prQS <- prQnu * cov/100 ## define starting values for prQS

				}
				else { ## Bensmail et al.
					g0 <- 0.5 + (r - 1)/2	## in general g0 must be a multiple of 0.5 for the 
								## Inverse Wishart (IW) to lead to a proper prior
					G0 <- g0 * solve(prQS)  ## match hierarchical and non-hierarchical priors
				}
				
				if(r > 1) {
					par.sigma <- list(c = array(prQnu, dim = c(1, K)),
							C = array(prQS, dim = c(r, r, K)),
							logdetC = array(detprQS, dim = c(1, K)),
							g = g0, G = G0)
				}
				else { ## r == 1
					par.sigma <- list(c = array(prQnu, dim = c(1, K)),
							C = array(prQS, dim = c(1, K)),
							g = g0, G = G0)
				}
			}
			else { ## non-hierarchical prior
				if(r > 1) {
					par.sigma <- list(c = array(prQnu, dim = c(1, K)),
							C = array(prQS, dim = c(r, r, K)),
							logdetC = array(detprQS, dim = c(1, K)))
				}
				else { ## r == 1 
					## later distinguish between 'sigmauniform' and 'others' ##
					par.sigma <- list(c = array(prQnu, dim = c(1, K)),
							C = array(prQS, dim = c(1, K)))
				}
			}
			par <- list(mu = par.mu, sigma = par.sigma)	
		}

		## prior degrees of freedom student-t ##
		if(dist == "student" || dist == "studmult") {
			## default prior: independent hierarchical prior following Fernandéz and Steel (1999)
			df.type <- "inhier" 
			df.trans <- 1
			df.a0 <- 2
			df.b0 <- 2
			df.mean <- 10
			df.d <- (df.mean - df.trans) * (df.b0 - 1)	
			df <- list(type = df.type, trans = df.trans, a0 = df.a0, b0 = df.b0, d = df.d)
			par <- list(par, df = df)

		} 

		## prior weights ##
		if(K > 1) {
			e0 <- 4
			weight <- matrix(e0, nrow = 1, ncol = K)
		}
		else { ## K = 1
			weight <- matrix()
		}
		prior <- new("prior", weight = weight, par = par, type = type, hier = hier)
		return(prior)
		
}
	
setMethod("show", "prior", function(object) {
					cat("Object 'prior'\n")
					cat("	type		:", class(object), "\n")
					cat("	hier		:", getHier(object), "\n")
					cat("	type (prior)	:", getType(object), "\n")
					cat("	par		: List of ", length(names(getPar(object))), "\n")
					if(!all(is.na(object@weight))) {
						cat("	weight		:", paste(dim(object@weight), collapse = "x"),"\n")
					}
				}
)
## Getters ##
## Generic set in 'model' class ##
setMethod("getWeight", "prior", function(.Object) {
					return(.Object@weight)
				}
) 
## Generic set in 'model' class ##
setMethod("getPar", "prior", function(.Object) {
					return(.Object@par)
				}
)
##setGeneric("getType", function(.Object) standardGeneric("getType"))
setMethod("getType", "prior", function(.Object) {
					return(.Object@type)
				}
)
setGeneric("getHier", function(.Object) standardGeneric("getHier"))
setMethod("getHier", "prior", function(.Object){
					return(.Object@hier)
				}
)
## R usual setters ##
## Generic set in 'model' class ##
setReplaceMethod("setWeight", "prior", function(.Object, value) {
						.Object@weight <- value
						validObject(.Object)
						return(.Object)
					}
)
## Generic set in 'model' class ##
setReplaceMethod("setPar", "prior", function(.Object, value) {
						.Object@par <- value
						validObject(.Object)
						return(.Object)
					}
)
## Generic set in 'model' class ##
setReplaceMethod("setType", "prior", function(.Object, value) {
						.Object@type <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setHier<-", function(.Object, value) standardGeneric("setHier<-"))
setReplaceMethod("setHier", "prior", function(.Object, value) {
						.Object@hier <- value
						validObject(.Object)
						return(.Object)
					}
)

