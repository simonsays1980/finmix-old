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
# along with finmix.  If not, see <http://www.gnu.org/licenses/>.

.prior <- setClass("prior",
                   representation(weight 	= "matrix",
                                  par 	    = "list",
                                  type 	    = "character",
                                  hier 	    = "logical"),
                   validity = function(object) {
                       type.choices <- c("condconjugate", "independent")
                       if (!(object@type %in% type.choices)) {               
                           stop("Unknown prior 'type'. 'type' must be 'independent' 
                                or 'condconjugate'.")
                       }
                       ## else: OK
                       TRUE		
                   },
                   prototype(weight   = matrix(),
                             par      = list(),
                             type     = character(),
                             hier     = logical()
                             )
)

"prior" <- function(weight = matrix(), par = list(), 
                    type = "independent", hier = TRUE) 
{
    object <- new("prior", weight = weight, par = par, 
                  type = type, hier = hier)
    return(object)
}

"priordefine" <- function(data = data(), model = model(), 
                          coef.mat = NULL, varargin = NULL) 
{

    .valid.Data.Obj(data, model@dist) 
    validObject(model)
    validObject(varargin)
    if (model@dist == "cond.poisson") {
        .valid.Coef.mat(model, coef.mat) 
    }
    object <- .prior(hier = TRUE, type = "independent")
    generatePrior(object, data = data, model = model, 
                  varargin = varargin, coef.mat = coef.mat)
}

setMethod("generatePrior", "prior", 
          function(object, data, model, varargin, coef.mat) 
          {
              dist <- model@dist
              if (dist == "poisson") {
                  object <- .generatePriorPoisson(object, data, model, 
                                                  varargin = varargin)
              } else if (dist == "cond.poisson") {
                  object <- .generatePriorCondPoisson(object, data, model, 
                                                      varargin = varargin, 
                                                      coef.mat)
              } else if (dist == "binomial") {
                  object <- .generatePriorBinomial(object, model)
              } else if (dist == "exponential") {
                  object <- .generatePriorExponential(object, data)
              } else {
                  object <- .generatePriorNorstud(object, data, model, 
                                                  varargin = varargin)
                  if (dist == "student" || dist == "studmult") {
                      object <- .generateDFPrior(object, model)
                  }
              }
              .generatePriorWeight(object, model)
          }
)

setMethod("show", "prior", 
          function(object) {
              cat("Object 'prior'\n")
              cat("     class       :", class(object), "\n")
              cat("     hier        :", object@hier, "\n")
              cat("     type        :", object@type, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              if(!all(is.na(object@weight))) {
                  cat("     weight      :",
                      paste(dim(object@weight), collapse = "x"), "\n")
              }		
          }
)
## Getters ##
## Generic set in 'model' class ##
setMethod("getWeight", "prior", 
          function(object) 
          {
              return(object@weight)
          }
) 
## Generic set in 'model' class ##
setMethod("getPar", "prior", 
          function(object) 
          {
              return(object@par)
          }
)
##setGeneric("getType", function(object) standardGeneric("getType"))
setMethod("getType", "prior", 
          function(object) 
          {
              return(object@type)
          }
)
setGeneric("getHier", function(object) standardGeneric("getHier"))
setMethod("getHier", "prior", 
          function(object)
          {
              return(object@hier)
          }
)
## R usual setters ##
## Generic set in 'model' class ##
setReplaceMethod("setWeight", "prior", 
                 function(object, value) 
                 {
                     object@weight <- value
                     validObject(object)
                     return(object)
                 }
)
## Generic set in 'model' class ##
setReplaceMethod("setPar", "prior", 
                 function(object, value) 
                 {
                     object@par <- value
                     validObject(object)
                     return(object)
                 }
)
## Generic set in 'model' class ##
setReplaceMethod("setType", "prior", 
                 function(object, value) 
                 {
                     object@type <- value
                     validObject(object)
                     return(object)
                 }
)
setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))
setReplaceMethod("setHier", "prior", 
                 function(object, value) 
                 {
                     object@hier <- value
                     validObject(object)
                     return(object)
                 }
)

### Private functions
### These functions are not exported 
".generatePriorPoisson" <- function(object, data.obj, 
                                    model.obj, varargin) 
{
    K <- model.obj@K
    if (data.obj@bycolumn) {
        datam <- data.obj@y
    } else {
        datam <- t(data.obj@y)
    }
    if (is.null(varargin)) {
        object@hier <- TRUE ## default prior is hierarchical
        object@type <- "condconjugate"
    } else {
        object@hier <- varargin@hier
        object@type <- "condconjugate"
    }
    ## default prior based on matching moments ##

    ## choose level of overdispersion, depending on the 
    ## ratio overdispersion/mean^2 ##
    ## no idea: data-based choice 
    mean <- mean(datam, na.rm = TRUE) 
    over <- as.numeric(var(datam, na.rm = TRUE) - mean)
    if (over > 0) {
        a0 <- mean^2/over
    } else {
        a0 <- 10
    }
    if (object@hier) {
        g0 <- 0.5
        G0 <- mean * g0/a0
        b0 <- g0/G0
        par <- list(a = array(a0, dim = c(1, K)),
                    b = array(b0, dim = c(1, K)),
                    g = g0, G = G0)
    } else {
        b0 = a0/mean
        par <- list(a = array(a0, dim = c(1, K)),
                    b = array(b0, dim = c(1, K)))
    }
    object@par <- par
    return(object)
}

".generatePriorCondPoisson" <- function(object, data.obj, model.obj, 
                                        varargin, coef.mat)
{
    K       <- model.obj@K
    if (data.obj@bycolumn) {
        datam <- data.obj@y
    } else {
        datam <- t(data.obj@y)
    }
    if (is.null(varargin)) {
        object@hier <- TRUE ## default prior is hierarchical
        object@type <- "condconjugate"
    } else {
        object@hier <- varargin@hier
        object@type <- varargin@hier 
    }
    if (object@type == "trunc.norm") {
        s0      <- sd(datam)
        par     <- list(s = array(s0, dim = c(1,K)))
    } else if (object@type == "condconjugate") {
        mean    <- mean(datam, na.rm = TRUE)
        over    <- var(datam, na.rm = TRUE) - mean
        if (over > 0) {
            a0 <- mean^2/over
        } else {
            a0  <- 10
        }
        if (object@hier) {
            g0          <- .5
            G0          <- mean * g0/a0
            be          <- g0/G0
            par         <- list(a = array(a0, dim = c(1, K)),
                                b = array(be, dim = c(1, K)),
                                g = g0, G = G0)
            par$a       <- par$a / 2^(rowSums(coef.mat) - 1)
            par$coef    <- coef.mat
        } else {
            be          <- a0/mean
            par         <- list(a = array(a0, dim = c(1, K)),
                                b = array(be, dim = c(1, K)))
            par$a       <- par$a / 2^(rowSums(coef.mat) - 1)
            par$coef    <- coef.mat
        }
    }
    object@par <- par
    return(object)
}

".generatePriorBinomial" <- function(object, model) 
{
    K           <- model@K
    object@type <- "condconjugate"
    ## uniform prior ##
    a0 <- 1
    b0 <- 1
    object@par <- list(a = array(a0, dim = c(1, K)),
                       b = array(b0, dim = c(1, K)))
    return(object)
}

".generatePriorExponential" <- function(object, data) 
{
    if (data@bycolumn) {
        datam <- data@y
    } else {
        datam <- t(data@y)
    }
    ## prior following Wagner (2007) ##
    object@hier <- FALSE
    object@type <- "condconjugate"
    a0 <- 0.1
    be <- mean(datam) * a0
    object@par <- list(a = array(a0, dim = c(1, K)), 
                b = array(be, dim = c(1, K))) 
    return(object)
}

".generatePriorNorstud" <- function(object, data, 
                                    model, varargin) 
{
    r   <- data@r
    K   <- model@K
    ## check if varargin is non-empty and prior object ##
    ## set hierarchical or non-hierarchical prior ##
    if(is.null(varargin)) { 
        ## default prior: independent hierarchical prior ##
        object@hier <- TRUE
        object@type <- "independent"
    } else {
        object@hier <- varargin@hier
        object@type <- varargin@type
    }
    conjugate.prior <- object@type == "condconjugate"
    bensmail <- FALSE      
    rich.green <- FALSE
    if (conjugate.prior || !hier) {
        bensmail <- TRUE ## Prior following Bensmail et al. 
    } else {
        rich.green <- TRUE ## Prior following Richardson and Green for r = 1
                           ##                 Stephens             for r = 2 only 
    }
    if (rich.green) {
        ## row vectors: dimension 1 x r
        max <- apply(datam, 2, max, na.rm = TRUE)
        min <- apply(datam, 2, min, na.rm = TRUE)
        mean <- (max + min) * .5
        cov <- diag((max - min)^2)
    } else {
        ## row vectors: dimension 1 x r
        mean <- apply(datam, 2, mean, na.rm = TRUE)
        cov <- ifelse(r > 1, cov(datam), var(datam))
    }		
    b0 <- mean
    if (conjugate.prior) {
        B0sc <- 1 ## info contained in a standard conjugate (sc) prior (equal to N0)
    } else {
        B0inv <- solve(cov) ## info contained in a non-conjugate prior,
						    ## i.e. either by Richardson Green or by Benmail et al. 
    }
    if (!conjugate.prior) {
        if (r > 1) {
            par.mu <- list(b = array(t(b0), dim = c(r, K)),
                           Binv = array(B0inv, dim = c(r, r, K)))
        } else { ## r = 1
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
    ## where:	
    ##      prQnu:	    degrees of freedom for Wishart and shape for Gamma
    ## 		prQS :  	shape for Q
    ## 	
    ## Select Q0 the prior mean of Q. 
    ## Determine prQS from prQS = Q0 * (prQnu - (r + 1)/2). This matches Q0 to the mean 
    ## of the Inverse Gamma or the Inverse Wishart distribution and to the mode of Q {-1}
    ## i.e. the Gamma and Wishart distribution respectively. 
    ## Further, variance shrinkage towards the ratio prQS/dfQpr, where dfQpr bounds the 
    ## ratio of the variances.

    dfQpr <- 2.5             ## this bounds the ratio of variances to 10 for r = 1
    prQnu <- dfQpr + (r-1)/2

    if (K == 1) {
        phi <- 1 ## c0 heterogeneity 
    } else {
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
    if (r > 1) {
        detprQS <- log(det(matrix(prQS)))
    } 
    if (hier) {
        if (rich.green) {
            if (r == 1) {
                g0 <- 0.2  ## Richardson and Green. Sampling from Gamma allows
                ## arbitrary g0: 
                ## WARNING: seems to cause problems in bayesf
                prQnu <- 2 ## Note that prQnu standard is changed here
            } else if (r == 2){
                g0 <- 0.3  ## Stephens	
                ## WARNING:  seems to cause problems in bayesf
                prQnu <- 3 ## prQnu is changed also in relation from standard
            } else { ## r > 2
                g0 <- 0.5 + (r - 1)/2
            }
            g0 <- 0.5 + (r - 1)/2
            G0 <- 100 * g0/prQnu * solve(cov) ## Stephens
            prQS <- prQnu * cov/100 ## define starting values for prQS           
        } else { ## Bensmail et al.
            g0 <- 0.5 + (r - 1)/2	## in general g0 must be a multiple of 0.5 for the 
            ## Inverse Wishart (IW) to lead to a proper prior
            G0 <- g0 * solve(prQS)  ## match hierarchical and non-hierarchical priors
        }
        if (r > 1) {
            par.sigma <- list(c = array(prQnu, dim = c(1, K)),
                              C = array(prQS, dim = c(r, r, K)),
                              logdetC = array(detprQS, dim = c(1, K)),
                              g = g0, G = G0)
        } else { ## r == 1
            par.sigma <- list(c = array(prQnu, dim = c(1, K)),
                              C = array(prQS, dim = c(1, K)),
                              g = g0, G = G0)
        }
    } else { ## non-hierarchical prior
        if (r > 1) {
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
    object@par <- list(mu = par.mu, sigma = par.sigma)	
    return(object)
}

".generateDfPrior" <- function(object) 
{
    ## default prior: independent hierarchical prior following Fernandéz and Steel (1999)
    df.type <- "inhier" 
    df.trans <- 1
    df.a0 <- 2
    df.b0 <- 2
    df.mean <- 10
    df.d <- (df.mean - df.trans) * (df.b0 - 1)	
    df <- list(type = df.type, trans = df.trans, a0 = df.a0, b0 = df.b0, d = df.d)
    object@par$df <- df        
    return(object)
}

".generatePriorWeight" <- function(object, model) 
{
    K   <- model@K
    if(K > 1) {
        e0 <- 4
        weight <- matrix(e0, nrow = 1, ncol = K)
    }
    else { ## K = 1
        weight <- matrix()
    }
    return(object)
}

### Validity 
### The coefficient matrix 'coef.mat' for 'cond.poisson'
### distributions with conditional prior must be a lower
### triangular matrix with ones on its diagonal.
### Further it must be of type 'matrix' or 'array' with
### dimension K x K.
".valid.Coef.mat" <- function(model.obj, coef.mat) {
    K <- model.obj@K
    if (is.null(coef.mat)) {
        stop("For a conditional Poisson mixture a coefficient matrix 
             'coef.mat' has to be provided.") 
    } else if (!is.null(coef.mat)) {
        if (!is.matrix(coef.mat) && !is.array(coef.mat)) {
            stop("Argument 'coef.mat' must be of type 'matrix' or 'array'.")
        } else if (nrow(coef.mat) != ncol(coef.mat)) {
            stop("Argument 'coef.mat' must be a quadratic 'matrix' or 'array'.")
        } else if (nrow(coef.mat) != K || ncol(coef.mat) != K) {
            stop("Dimension of argument 'coef.mat' must correspond to number 
                 of components 'K' in 'model'.\n")
        } else if (!(all(diag(coef.mat) == 1))) {
            stop("Coefficients on the diagonal of 'coef.mat' must be equal 
                 to one.\n")
        }
    }
}

### 'data' objects must have a slot @y with not all values 
### NA. Further the dimensions @N and @r must be conform
### with the slot @y given the slot @bycolumn.
### Furthermore, only Student-t and Normal models are
### intended to handle multivariate data.
".valid.Data.Obj" <- function(data.obj, model.dist)
{
    validObject(data.obj)
    if (all(is.na(data.obj@y))) {
        stop("Argument 'data': slot 'y' is empty.")
    } else {
        if (data.obj@bycolumn) {
            N   <- nrow(data.obj@y)
            r   <- ncol(data.obj@y)
            if (N != data.obj@N) {
                data.obj@N  <- N
            }
            if (r != data.obj@r) {
                data.obj@r  <- r
            }
        } else {
            N   <- ncol(data.obj@y)
            r   <- nrow(data.obj@y)
            if (N != data.obj@N) {
                data.obj@N  <- N
            }
            if (r != data.obj@r) {
                data.obj@r  <- r
            }
        }
        if (data.obj@r > 1) {
            if (!model.dist %in% c("normult", "studmult")) {
                stop("Model", model.dist, "can not handle multivariate
                     data")
            } 
        }
    }
}


