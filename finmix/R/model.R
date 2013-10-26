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
# along with finmix. If not, see <http://www.gnu.org/licenses/>.

.model <- setClass("model",
                   representation(dist        ="character",
                                  r           = "integer",
                                  K           = "integer",
                                  weight      = "matrix",
                                  par         = "list",
                                  indicmod    = "character",
                                  indicfix    = "logical",
                                  T           = "matrix"),
                   validity = function(object) {
                       .valid.Model(object)
                       ## else: OK ##
                       TRUE
                   },
                   prototype(dist     = character(),
                             r        = integer(),
                             K        = integer(),
                             weight   = matrix(),
                             par      = list(),
                             indicmod = character(),
                             indicfix = logical(),
                             T        = matrix()
                             )
)

## Constructor for class 'model' ##
"model" <- function(dist = "poisson", r, K, 
                    weight = matrix(), par = list(), 
                    indicmod = "multinomial", 
           			indicfix = FALSE, T = matrix()) 
{
    if (missing(K)) {
        K <- .check.K.Model(weight)
    } else {
        K <- as.integer(K)
    }
    if (missing(r)) {
        r <- .check.r.Model(dist) 
    } else {
        r <- as.integer(r)
    }
    if (missing(weight) && K > 1) {
        weight <- .check.weight.Model(K)
    } else {
        weight <- as.matrix(weight)
    }
    if (!missing(T)) {
        T <- .check.T.Model(T)
    }
    .model(dist = dist, r = r, K = K, weight = weight, 
           par = par, indicmod = indicmod, 
           indicfix = indicfix, T = T)
}

setMethod("hasWeight", "model",
          function(object, verbose = FALSE) 
          {
              if (!all(is.na(object@weight))) {
                  if (ncol(object@weight) == object@K) {                      
                      return(TRUE)
                  } else {
                      if (verbose) {
                          stop(paste("Wrong dimension of ",
                                     "slot 'weight' of ",
                                     "'model' object." ,
                                     "Weights must be of ",
                                     "dimension 1 x K.",
                                     sep = ""))
                      } else {
                          return(FALSE)
                      }
                  }
              } else {
                  if (verbose) {
                      stop(paste("Slot 'weight' of 'model' ",
                                 "object is empty.",
                                 sep = ""))                      
                  } else {                      
                      return(FALSE)
                  }
              }
          }
)

setMethod("hasT", "model",
          function(object, verbose = FALSE) 
          {
              if (!all(is.na(object@T))) {
                  return(TRUE)
              } else {
                  if (verbose) {
                      stop(paste("Slot 'T' of 'model' ",
                                 "object is empty.", 
                                 sep = ""))
                  } else {
                      return(FALSE)
                  }
              }
          }
)

setMethod("hasPar", "model",
          function(object, verbose = FALSE) 
          {
              .haspar.Model(object, verbose)
          }
)

### ----------------------------------------------------------------------
### Simulate method
### @description    Simulates values for a specified model in an 'model'
###                 object.
### @par    model       an S4 'model' object; with specified parameters 
### @par    N           an R 'integer' value specifying the number of
###                     values to be simulated
### @par    varargin    an S4 'fdata' object; with specified variable 
###                     dimension @r and repetitions @T
### @return         an S4 object of class 'fdata' holding the simulated
### @see    ?simulate
### @author Lars Simon Zehnder
### ----------------------------------------------------------------------
setMethod("simulate", "model", 
          function(model, N = 100, varargin, seed = 0)
          {
## TODO: CHeck model for parameters. Check varargin for dimension. Check
##      model anf varargin for consistency.
              if (!missing(seed)) {
                  set.seed(seed)
              } ## Implemented maybe finmixOptions with a state variable seed
              if (!hasWeight(model)) {
                  model@weight <- matrix(1/model@K, nrow = 1, ncol = model@K)
              }
              ## Start simulating the allocations
              S <- .simulate.indicators.Model(model, N) 
              if (missing(varargin)) {
                  varargin <- fdata(r = 1, T = matrix(1, nrow = N), S = S)
              } else {
                  varargin@S <- S
              }
              .simulate.data.Model(model, N, varargin)
          }
)

## plot ##
setMethod("plot", "model", 
          function(x, y, dev = TRUE, ...) 
          {
              dist    <- x@dist
              if(dist == "normal") {
                  .plot.Normal.Model(x, dev, ...)
              } else if (dist == "normult") {
                  .plot.Normult.Model(x, dev, ...)
              } else if (dist == "exponential") {
                  .plot.Exponential.Model(x, dev, ...)
              } else if (dist == "student") {
                  .plot.Student.Model(x, dev, ...)
              } else if (dist == "studmult") {
                  .plot.Studmult.Model(x, dev, ...)
              }	else if (dist %in% c("poisson", "cond.poisson")) {
                  .plot.Poisson.Model(x, dev, ...)
              } else if (dist == "binomial") {
                  if(length(x@T) != 1) {
                      stop("Plotting a binomial distribution with varying 
                           repetitions in slot 'T' is not possible.")
                  }
                  .plot.Binomial.Model(x, dev, ...)
              }
          }
)

setMethod("plotPointProc", signature(x      = "model",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              hasPar(x, verbose = TRUE)
              hasWeight(x, verbose = TRUE)
              if (x@dist == "poisson") {
                  .plotpointproc.Poisson(x, dev)
              }
          }
)

## Marginal Mixture ##
setMethod("mixturemar", "model", 
          function(object, J) 
          {
              .mixturemar.Model(object, J)
          }
)
	
## Show ##
setMethod("show", "model", 
          function(object) 
          {
              cat("Object 'model'\n")
              cat("     class       :", class(object), "\n")
              cat("     dist        :", object@dist, "\n")
              cat("     r           :", object@r, "\n")
              cat("     K           :", object@K, "\n")
              if (hasPar(object)) {
                 cat("     par         : List of",
                      length(object@par), "\n")
              }
              if (!object@indicfix) {
                  cat("     weight      :", 
                      paste(dim(object@weight), collapse = "x"), 
                      "\n")
              }
              cat("     indicmod    :", object@indicmod, "\n")
              cat("     indicfix    :", object@indicfix, "\n")
              if (object@dist == "binomial" && !all(is.na(object@T))) {
                  cat("     T           :",
                      paste(dim(object@T), collapse = "x"), "\n")
              }
          }
)

## Getters ##
setMethod("getDist", "model", 
          function(object) 
          {
              return(object@dist)
          }
)

setMethod("getR", "model", 
          function(object) 
          {
              return(object@r)
          }
)

setMethod("getK", "model", 
          function(object) 
          {
              return(object@K)
          }
)

setMethod("getWeight", "model", 
          function(object) 
          {
              return(object@weight)
          }
)

setMethod("getPar", "model", 
          function(object) 
          {
              return(object@par)
          }
)

setMethod("getIndicmod", "model", 
          function(object) 
          {
              return(object@indicmod)					
          }
)

setMethod("getIndicfix", "model", 
          function(object) 
          {
              return(object@indicfix)
          }
)

setMethod("getT", "model", 
          function(object) 
          {
              return(object@T)
          }
)

## Setters ##
setReplaceMethod("setDist", "model", 
                 function(object, value) 
                 {
                     object@dist <- value
                     .valid.dist.Model(object)
                     return(object)
                 }
)

setReplaceMethod("setR", "model", 
                 function(object, value)
                 {
                     object@r <- as.integer(value)
                     validObject(object)
                     return(object)
                 }
)

setReplaceMethod("setK", "model", 
                 function(object, value) 
                 {
                     object@K           <- as.integer(value)
                     .init.valid.K.Model(object)
                     if (object@K > 1) {
                         object@weight  <- .check.weight.Model(object@K)                   
                     } else {
                         weight                 <- matrix()
                         storage.mode(weight)   <- "numeric"
                         object@weight          <- weight
                     }
                     return(object)
                 }
)

setReplaceMethod("setWeight", "model", 
                 function(object, value) 
                 {
                     object@weight  <- as.matrix(value)
                     object@K       <- ncol(object@weight)
                     .init.valid.weight.Model(object)
                     return(object)
                 }
)

setReplaceMethod("setPar", "model", 
                 function(object, value) 
                 {
                     object@par <- value
                     .init.valid.par.Model(object)                     
                     return(object)
                 }
)

setReplaceMethod("setIndicmod", "model", 
                 function(object, value) 
                 {
                     object@indicmod <- value
                     return(object)
                 }
)

setReplaceMethod("setIndicfix", "model", 
                 function(object, value) 
                 {
                     object@indicfix <- value
                     return(object)
                 }
)

setReplaceMethod("setT", "model",
                 function(object, value) 
                 {                     
                     object@T <- matrix(value)                    
                     .valid.T.Model(object)
                     return(object)
                 }
)

### Private functions
### These functions are not exported

### Checking.
### Checking is used for in the constructor. 
### Arguments for the slots are checked for validity and
### if missing are given by default values. Altogether the
### constructor tries to construct a fully specified model
### object with consistent slots.

### Check K: If weights are provided by the user, the number
### of components is set to the number of columns of the weights.
### If argument 'weight' is missing from the call, the number of
### components is assumed to be one. 
".check.K.Model" <- function(weight)
{
    if (!all(is.na(weight))) {
        return(NCOL(weight))
    } else {
        return(as.integer(1))
    }
}

### Check r: The dimension of the model is determined in regard to 
### the defined distribution in argument 'dist' (if missing the 
### default is 'poisson'). For univariate distributions it is set
### to one and for multivariate distribution as a default to two.
".check.r.Model" <- function(dist) 
{
    univ    <- .get.univ.Model()
    multiv  <- .get.multiv.Model()
    if (dist %in% univ) {
        return(as.integer(1))
    } else if (dist %in% multiv) {
        return(as.integer(2))
    } else {
        stop(paste("Unknown distribution in slot ",
                   "'dist' of 'model' object.", 
                   sep = ""))
    }
}

### Check weight: If argument 'weight' is missing from the call 
### equally balanced weights are given as a default. 
".check.weight.Model" <- function(K) 
{
    weight <- matrix(1/K, nrow = 1, ncol = K)
    return(weight)
} 

### Check T: If repetitions are given they are checked in regard
### to validity. In case of non-numeric objects an error is thrown.
### In case of objects of type 'numeric' it is implicitly converted 
### to type 'integer'.
".check.T.Model" <- function(T)
{
    if (!all(is.na(T))) {
        if (!is.numeric(T)) {
            stop(paste("Wrong specification of slot 'T' in ",
                       "'model' object. Repetitions must be of ",
                       "type 'integer'.", sep = ""))
        } else {
            storage.mode(T) <- "integer"
            return(T)
        }
    }
}

### Marginal model
".mixturemar.Model" <- function(obj, J)
{
    if (object@dist == "normult") {
        .mixturemar.normult.Model(obj, J)
    } else if (object@dist == "studmult") {
        .mixturemar.studmult.Model(obj, J)
    } else {
        stop("A marginal distribution can only be obtained from 
             multivariate distributions.")
    }
}

".mixturemar.normult.Model" <- function(obj, J)
{
    dist <- ifelse(length(J) == 1, "normal", "normult")
    r             <- length(J)
    K             <- obj@K
    weight        <- obj@weight
    mu            <- obj@par$mu[J, ]
    sigma         <- obj@par$sigma[J, J, ]
    par           <- list(mu = mu, sigma = sigma)
    indicmod      <- "multinomial"
    indicfix      <- TRUE
    margin.model  <- .model(dist = dist, r = r, K = K, 
                            weight = weight, par = par, 
                            indicmod = indicmod, 
                            indicfix = indicfix)
    validObject(margin.model)
    return(margin.model)
}

".mixturemar.studmult.Model" <- function(obj, J)
{
    dist <- ifelse(length(J) == 1, "student", "studmult")
    r             <- length(J)
    K             <- obj@K
    weight        <- obj@weight
    mu            <- obj@par$mu[J, ]
    sigma         <- obj@par$sigma[J, J, ] 
    df            <- obj@par$df
    par           <- list(mu = mu, sigma = sigma, df = df)
    indicmod      <- "multinomial"
    indicfix      <- TRUE
    margin.model  <- .model(dist = dist, r = r, K = K, 
                            weight = weight, par = par, 
                            indicmod = indicmod,
                            indicfix = indicfix)
    validObject(margin.model)
    return(margin.model)
}

### ==============================================================
### Simulate
### --------------------------------------------------------------

### --------------------------------------------------------------
### .simulate.indicators.Model
### @description    Simulates the indicators.
### @par    obj an S4 object of class 'model' 
### @par    N   an R 'integer' object
### @return         an R 'matrix' object with N simulated indi-
###                 cators.
### @details        indicators are simulated via the slot @weight
###                 the 'model' object
### @see    ?simulate
### @author Lars Simon Zehnder
### --------------------------------------------------------------

### TODO: Implement C++ function.
".simulate.indicators.Model" <- function(obj, N) 
{
    K <- obj@K
    if (K == 1) {
        S <- matrix(as.integer(1), nrow = N, ncol = K)
    } else {
        ## if (model@indicmod = "") -> "Multinomial"
        ## if Markov else
        if (obj@indicmod == "multinomial") {
            rnd       <- runif(N)
            rnd       <- matrix(rnd, nrow = N, ncol = 2)
            weightm   <- matrix(obj@weight, nrow = N, ncol = K, 
                                byrow = TRUE)
            S         <- apply((t(apply(weightm, 1, cumsum)) < rnd), 1, sum) + 1
            S         <- matrix(S, nrow = N)                      
        }
    }
    return(S)
}

### --------------------------------------------------------------------
### .simulate.data.Model
### @description    Simulates the simulation functions for a specific model.
### @par    obj         an S4 'model' object
### @par    N           an R 'integer' object; number of simulated values
### @par    fdata.obj   an S4 'fdata' object
### @return         an S4 object of class 'fdata' with simulated values
### @see    ?fdata, ?simulate
### @author Lars Simon Zehnder
### ---------------------------------------------------------------------
".simulate.data.Model" <- function(obj, N, fdata.obj)
{
    dist <- obj@dist
    if (dist == "poisson") {
        .simulate.data.poisson.Model(obj, N, fdata.obj)
    } else if (dist == "binomial") {
        .simulate.data.binomial.Model(obj, N, fdata.obj)
    }
}

### ---------------------------------------------------------------------
### .simulate.data.poisson.Model
### @description    Simulates values from a Poisson mixture using pre-
###                 specified model and indicators
### @par    obj         an S4 object of class 'model'
### @par    N           an R 'integer' object; number of simulated values
### @par    fdata.obj   an S4 object of class 'fdata'
### @return         an S4 object of class 'fdata' with simulated values
### @see    ?simulate, model:::.simulate.data.Model, ?rpois
### @author Lars Simon Zehnder
### ---------------------------------------------------------------------
".simulate.data.poisson.Model" <- function(obj, N, fdata.obj)
{
    fdata.obj@type  <- "discrete"
    fdata.obj@sim   <- TRUE
    fdata.obj@y         <- matrix(rpois(N, obj@par$lambda[fdata.obj@S]))
    return(fdata.obj)
   
}

### ---------------------------------------------------------------------
### .simulate.data.binomial.Model
### @description    Simulates values from a Binomial mixture using pre-
###                 specified model and indicators
### @par    obj         an S4 object of class 'model'
### @par    N           an R 'integer' object; number of simulated values
### @par    fdata.obj   an S4 object of class 'fdata'
### @return         an S4 object of class 'fdata' with simulated values
### @see    ?simulate, model:::.simulate.data.Model, ?rbinom
### @author Lars Simon Zehnder
### ---------------------------------------------------------------------
".simulate.data.binomial.Model" <- function(obj, N, fdata.obj)
{
    if (!hasT(fdata.obj)) {
        fdata.obj@T <- as.matrix(1)
    }
    fdata.obj@type  <- "discrete"
    fdata.obj@sim   <- TRUE
    fdata.obj@y     <- matrix(rbinom(N, fdata.obj@T, obj@par$p[fdata.obj@S]))    
    return(fdata.obj)
}

### Plotting
### Plot Poisson models: Poisson models are discrete
### models and a barplot is used. 
### The range for the x-axis is determined via the 
### quantiles of the largest and smallest Poisson model
### in the mixture. 
".plot.Poisson.Model" <- function(model.obj, dev, ...)
{

    if (.check.grDevice() && dev) {
        dev.new(title = "Model plot")
    }
    lambda      <- model.obj@par$lambda
    weight      <- model.obj@weight
    xlim.up     <- qpois(.9999, lambda = max(lambda))
    xlim.low    <- qpois(.0001, lambda = min(lambda))
    x.grid      <- seq(xlim.low, xlim.up, by = 1)
    y.grid      <- sapply(x.grid, dpois, lambda = lambda)
    y.grid      <- weight %*% y.grid
    main.title  <- paste("Poisson Mixture K = ", 
                         model.obj@K, sep="")
    label.grid  <- axisTicks(c(xlim.low, xlim.up), log = FALSE,
                             nint = 10)
    bp          <- barplot(y.grid, main = main.title, axes = F,  
                           col = "gray65", border = "gray65", ...)
    axis(side = 2, cex = .7, cex.axis = .7)
    axis(side = 1, tick = FALSE, at = bp[which(x.grid %in% label.grid)], 
         labels = label.grid, cex.axis = .7)
    mtext(side = 1, "x", cex = .7, cex.axis = .7, line = 3)
    mtext(side = 2, "P(x)", cex = .7, cex.axis = .7, line = 3)
}

### Plot Binomial models: Binomial models are discrete 
### models and line model is used. 
### The grid for the x-axis is determined by taking 
### the 
".plot.Binomial.Model" <- function(model.obj, dev, ...)
{
    if (.check.grDevice() && dev) {
        dev.new(title = "Model plot")
    }
    n           <- model.obj@T
    p           <- model.obj@par$p
    weight      <- model.obj@weight
    xlim        <- max(n, na.rm = TRUE)
    x.grid      <- seq(0, xlim, by = 1)
    y.grid      <- sapply(x.grid, dbinom, size = n, p = p)
    y.grid      <- weight %*% y.grid
    main.title  <- paste("Binomial Mixture K = ", 
                         model.obj@K, sep = "")
    plot(x.grid, y.grid, main = main.title, type = "h", 
         xlab = "x", ylab = "P(x)", ...)
    points(x.grid, y.grid, pch = 20)
}

".plot.Exponential.Model" <- function(model.obj, dev, ...)
{
    if (.check.grDevice() && dev) {
        dev.new(title = "Model plot") 
    }
    lambda      <- model.obj@par$lambda
    weight      <- model.obj@weight
    min.lambda  <- min(lambda, na.rm = TRUE)
    xlim        <- qexp(.9999, rate = min.lambda)
    x.grid      <- seq(0, ceiling(xlim), length = 
                       as.integer(100 * lambda^(-2)))
    y.grid      <- sapply(x.grid, dexp, rate = lambda)
    y.grid      <- weight %*% y.grid
    main.title  <- paste("Exponential Mixture K = ",
                         model.obj@K, sep = "")
    plot(x.grid, y.grid, main = main.title, type = "l", 
         xlab = "x", ylab = "P(x)", ...)
}

".plot.Student.Model" <- function(model.obj, dev, ...)
{
    if (.check.grDevice() && dev) {
        dev.new(title = "Model plot") 
    }
    mu          <- model.obj@par$mu
    sigma       <- model.obj@par$sigma
    df          <- model.obj@par$df
    weight      <- model.obj@weight
    max.mu      <- max(mu, na.rm = TRUE)
    max.sigma   <- max(sigma, na.rm = TRUE)
    min.df      <- min(df, na.rm = TRUE)
    xlim        <- max.mu + max.sigma * qt(.9999, min.df)
    x.grid      <- seq(-xlim, xlim, length = 1000) + max.mu
    y.grid      <- sapply(x.grid, "-", mu)
    y.grid      <- apply(y.grid, 2, "/", sigma)
    y.grid      <- apply(y.grid, 2, dt, df = df)
    y.grid      <- apply(y.grid, 2, "/", sqrt(sigma))
    y.grid      <- t(weight %*% y.grid)
    main.title  <- paste("Student-t Mixture K = ", 
                         model.obj@K, sep="")
    plot(x.grid, y.grid, main = main.title, type = "l", 
         xlab = "x", ylab = "P(x)", ...)
}

".plot.Normal.Model" <- function(model.obj, dev, ...)
{
    if (.check.grDevice() && dev) {
        dev.new(title = "Model Plot")
    }
    mu          <- model.obj@par$mu
    sigma       <- model.obj@par$sigma
    weight      <- model.obj@weight
    max.mu      <- max(mu, na.rm = TRUE)
    max.sigma   <- max(mu, na.rm = TRUE)
    xlim        <- qnorm(.9999, mean = max.mu, 
                         sd = max.sigma)
    x.grid      <- seq(-xlim, xlim, length = 1000) + max.mu
    y.grid      <- sapply(x.grid, dnorm, mean = mu, 
                          sd = sigma)
    y.grid      <- weight %*% y.grid
    main.title  <- paste("Normal Mixture K = ",
                         model.obj@K, sep = "")
    plot(x.grid, y.grid, main = main.title, type = "l",
         xlab = "x", ylab = "P(x)", ...)
}

".plot.Normult.Model" <- function(model.obj, dev, ...)
{
    K               <- model.obj@K
    r               <- model.obj@r
    if (r == 2) {
        if (.check.gr.Device() && dev) {
            dev.new(title = "Model: Perspective plot")
        }
        xyz.grid    <- .generate.Grid.Normal(model.obj)
        main.title = paste("Multivariate Normal Mixture K = ",
                           K, sep = "")
        persp(xyz.grid$x, xyz.grid$y, xyz.grid$z, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, expand = 0.5,
              lphi = 180, ltheta = 90, r = 40, d = 0.1, 
              ticktype = "detailed", zlab = "P(x)", xlab = "r = 1",
              ylab = "r = 2", cex = 0.7, cex.lab = 0.7, cex.axis = 0.7)
    } else if (r > 2 && r < 6) {
        if (.check.grDevice() && dev) {
            dev.new(title = "Model: Contour plots") 
        }
        if (r == 3) {
            par(mfrow = c(1, r), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
        } else if (r == 4) {
            par(mfrow = c(2, 3), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
        } else {
            par(mfrow = c(2, 5), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
        }
        for (i in seq(1, r - 1)) {
            for (j in seq(1, r)) {
                marmodel    <- mixturemar(model.obj, J = c(i, j))
                xyz.grid    <- .generate.Grid.Normal(marmodel)
                contour(xyz.grid$x, xyz.grid$y, xyz.grid$z, 
                        col = "gray47", cex = 0.7, cex.axis = 0.7,
                        xlab = paste("r = ", i, sep = ""), 
                        ylab = paste("r = ", j, sep = ""))
            }
        }
    } else {
        stop("Method 'plot' for 'model' objects is not implemented for
             model dimensions of r > 5.")
    }
}

".plot.Normult.Model" <- function(model.obj, dev, ...)
{
    K               <- model.obj@K
    r               <- model.obj@r
    if (r == 2) {
        if (.check.gr.Device() && dev) {
            dev.new(title = "Model: Perspective plot")
        }
        xyz.grid    <- .generate.Grid.Student(model.obj)
        main.title = paste("Multivariate Student-t Mixture K = ",
                           K, sep = "")
        persp(xyz.grid$x, xyz.grid$y, xyz.grid$z, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, expand = 0.5,
              lphi = 180, ltheta = 90, r = 40, d = 0.1, 
              ticktype = "detailed", zlab = "P(x)", xlab = "r = 1",
              ylab = "r = 2", cex = 0.7, cex.lab = 0.7, cex.axis = 0.7)
    } else if (r > 2 && r < 6) {
        if (.check.grDevice() && dev) {
            dev.new(title = "Model: Contour plots") 
        }
        if (r == 3) {
            par(mfrow = c(1, r), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
        } else if (r == 4) {
            par(mfrow = c(2, 3), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
        } else {
            par(mfrow = c(2, 5), mar = c(2, 2, 2, 2),
                oma = c(4, 5, 1, 5))
        }
        for (i in seq(1, r - 1)) {
            for (j in seq(1, r)) {
                marmodel    <- mixturemar(model.obj, J = c(i, j))
                xyz.grid    <- .generate.Grid.Student(marmodel)
                contour(xyz.grid$x, xyz.grid$y, xyz.grid$z, 
                        col = "gray47", cex = 0.7, cex.axis = 0.7,
                        xlab = paste("r = ", i, sep = ""), 
                        ylab = paste("r = ", j, sep = ""))
            }
        }
    } else {
        stop("Method 'plot' for 'model' objects is not implemented for
             model dimensions of r > 5.")
    }
}

".generate.Grid.Normal" <- function(model.obj)
{
    mu              <- model.obj@par$mu
    sigma           <- model.obj@par$sigma
    weight          <- model.obj@weight
    func <- function(s, t) 
    {
        value <- 0
        for (k in seq(1, K)) {
            value <- value + weight[k] * 
            dmvnorm(cbind(s, t), mean = mu[, k], 
                    sigma = sigma[,, k])
        }
    }    
    mu.norm         <- apply(mu, 2, function(x) sqrt(sum(x^2)))
    max.mu.index    <- tail(sort(mu.norm, index = TRUE)$ix, 1)
    max.mu          <- mu[, max.mu.index]
    sigma.det       <- apply(sigma, 3, det)
    max.sigma.index <- tail(sort(sigma.det, index = TRUE)$ix, 1)
    max.sigma       <- sigma[,, max.sigma.index]
    xylim           <- qmvnorm(.9999, mean = max.mu, 
                               sigma = max.sigma)$quantile
    x.grid          <- seq(-xylim, xylim, length = 100)
    xy.grid         <- cbind(x.grid, x.grid) 
    xy.grid         <- t(apply(xy.grid, 1, "+", max.mu))
    z.grid          <- outer(xy.grid[, 1], xy.grid[, 2], func)
    grid.list       <- list(x = xy.grid[, 1], y = y.grid[, 2], 
                            z = z.grid)
    return(grid.list)
}

".generate.Grid.Student" <- function(model.obj)
{
    mu              <- model.obj@par$mu
    sigma           <- model.obj@par$sigma
    df              <- model.obj@par$df
    weight          <- model.obj@weight
    func <- function(s, t) 
    {
        value <- 0
        for (k in seq(1, K)) {
            value <- value + weight[k] * 
            dmvt(cbind(s, t), delta = mu[, k], 
                 sigma = sigma[,, k], df = df[k])
        }
    }    
    mu.norm         <- apply(mu, 2, function(x) sqrt(sum(x^2)))
    max.mu.index    <- tail(sort(mu.norm, index = TRUE)$ix, 1)
    max.mu          <- mu[, max.mu.index]
    sigma.det       <- apply(sigma, 3, det)
    max.sigma.index <- tail(sort(sigma.det, index = TRUE)$ix, 1)
    max.sigma       <- sigma[,, max.sigma.index]
    min.df          <- min(df, na.rm = TRUE)
    xylim           <- qmvt(.9999, delta = max.mu, 
                            sigma = max.sigma, df = min.df)$quantile
    x.grid          <- seq(-xylim, xylim, length = 100)
    xy.grid         <- cbind(x.grid, x.grid) 
    xy.grid         <- t(apply(xy.grid, 1, "+", max.mu))
    z.grid          <- outer(xy.grid[, 1], xy.grid[, 2], func)
    grid.list       <- list(x = xy.grid[, 1], y = y.grid[, 2], 
                            z = z.grid)
    return(grid.list)
}

### plotPointProc
".plotpointproc.Poisson" <- function(x, dev) 
{
    K   <- x@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Point Process Representation")
    }
    if (min(x@par$lambda) < 1) {
        lambda <- log(x@par$lambda)
    } else {
        lambda <- x@par$lambda
    } 
    y.grid      <- rep(0, K)
    size.grid   <- as.vector(x@weight * 4)
    col.grid    <- gray.colors(K, start = 0.2,
                               end = 0.5)
    plot(lambda, y.grid, pch = 20, col = col.grid, 
         cex = size.grid, cex.lab = .7, cex.axis = .7, 
         main = "", ylab = "", xlab = "")
    mtext(side = 1, bquote(lambda), cex = .7, cex.lab = .7, 
          line = 3)
    legend.names    <- list("", K)
    for (k in seq(1, K)) {
        legend.names[[k]]   <- bquote(lambda[.(k)])
    }
    legend("topright", legend = do.call(expression, legend.names), 
           col = col.grid, fill = col.grid)
}

### Has
### Checks if a 'model' object has specified parameters.
".haspar.Model" <- function(obj, verbose) 
{
    if (length(obj@par) > 0) {
        dist <- obj@dist
        if (dist %in% c("poisson", "cond.poisson")) {
            .haspar.poisson.Model(obj, verbose)
        } else if (dist == "binomial") {
            .haspar.binomial.Model(obj, verbose)
        }
    } else {
        if (verbose) {
            stop(paste("Slot 'par' of 'model' object is ",
                       "empty.", sep = ""))
        } else {
            return(FALSE)
        }
    }
}

### -----------------------------------------------------------------
### .haspar.poisson.Mode
### @description    Checks if a Poisson model has fully specified 
###                 parameters. If verbose is set to TRUE an error 
###                 is thrown.
### @par    obj     an S4 object of class 'model'
### @par    verbose an object of class 'logical'
### @return         either TRUE or FALSE if parameters are fully 
###                 specified or not. In case verbose == FALSE an 
###                 error is thrown.
### -----------------------------------------------------------------
".haspar.poisson.Model" <- function(obj, verbose) 
{
    if ("lambda" %in% names(obj@par)) {
        if (length(obj@par$lambda) != obj@K) {
            if (verbose) {
                stop(paste("Wrong specification of slot @par of ",
                           "'model' object. Number of Poisson ",
                           "parameters in @par$lambda must match ",
                           "number of components in slot @K.", 
                           sep = ""))
            } else {
                return(FALSE)
            }
        } else {
            return(TRUE)
        }
    } else {
        if (verbose) 
        {
            stop(paste("Wrong specification of slot @par of ",
                       "'model' object. Poisson parameters must be ",
                       "named 'lambda'.", sep = ""))
        } else {
            return(FALSE)
        }
    }
}
### -------------------------------------------------------------------
### .haspar.binomial.Model
### @description    Checks if a Binomial model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @par    obj     an S4 object of class 'model'
### @par    verbose an object of class 'logical'
### @return         either TRUE or FALSE if parameters are fully 
###                 specified or not. In case verbose == TRUE an
###                 error is thrown.
### -------------------------------------------------------------------
".haspar.binomial.Model" <- function(obj, verbose)
{
    if ("p" %in% names(obj@par)) {
        if (length(obj@par$p) != obj@K) {
            if (verbose) {
                stop(paste("Wrong specification of slot @par of ",
                           "'model' object. Number of Binomial ",
                           "parameters in @par$p must match ",
                           "number of components in slot @K.",
                           sep = ""))
            } else {
                return(FALSE)
            }
        } else {
            return(TRUE)
        }
    } else {
        if (verbose) {
            stop(paste("Wrong specification of slot @par of ",
                       "'model' object. Binomial parameters must be ",
                       "named 'p'.", sep = ""))
        } else {
            return(TRUE)
        }
    }
}
### Validity
### Validity checking of model objects is implemented
### in two versions: an initializing version relying partly
### on warnings and amore restrictive version relying exclusively
### on errors. 
### The less restrictive validity check is used in setters and
### and the fully restrictive version in the constructor and later
### usage of model object (e.g. see 'mcmcstart()').
".init.valid.Model" <- function(obj) 
{
    .valid.dist.Model(obj)
    .init.valid.K.Model(obj) 
    .init..valid.r.Model(obj)
    .init.valid.par.Model(obj)
    .init.valid.weight.Model(obj)
    .valid.T.Model(obj)
}

".valid.Model" <- function(obj) 
{
    .valid.dist.Model(obj)
    .valid.K.Model(obj)
    .valid.r.Model(obj)
    .valid.par.Model(obj)
    .valid.weight.Model(obj)
    .valid.T.Model(obj)
}

### Valid dist: @dist must be a distribution which
### is implemented.
### Furthermore: @indicmod must be an indicator 
### model that is implemented.
".valid.dist.Model" <- function(obj)
{
    dists            <- c("normal", "normult", "exponential", 
                          "student", "studmult", "poisson", 
                          "cond.poisson", "binomial")
    indicmod.dists   <- c("multinomial")
    if (length(obj@dist) > 0) {
        if (!(obj@dist %in% dists)) {
            stop(paste("Unknown distribution in slot 'dist' ",
                          "of 'model' object.", sep = ""))            
        } else {
            if (!(obj@indicmod %in% indicmod.dists)) {
            stop(paste("Unknown indicator distribution in slot ", 
                          "'indicmod' of 'model' object."), sep = "")
            }
        }
    }
}

### Valid K: The number of components @K must be a positive 
### integer with no exclusion.
".init.valid.K.Model" <- function(obj) 
{
    if (obj@K < 1) {
        stop(paste("Wrong specification of slot 'K' of ",
                   "'model' object. Number of components ",
                   "must be a positive integer.", sep = ""))
    } else {
        if (!all(is.na(obj@weight))) {            
            if (obj@K != ncol(obj@weight)) {
                warning(paste("Dimension of slot 'weight' in ",
                              "'model' object does not match ",
                              "number of components in slot 'K'.",
                              sep = ""))
            }
        }
        .init.valid.par.Model(obj)            
     }

}

".valid.K.Model" <- function(obj) 
{
     if (obj@K < 1) {
        stop(paste("Wrong specification of slot 'K' of ",
                   "'model' object. Number of components ",
                   "must be a positive integer.", sep = ""))
    } else {
        if (!all(is.na(obj@weight))) {            
            if (obj@K != ncol(obj@weight)) {
                stop(paste("Dimension of slot 'weight' in ",
                           "'model' object does not match ",
                           "number of components in slot 'K'.",
                           sep = ""))
            }
        }
        .valid.par.Model(obj)            
     }
}
   
### Valid r: @r must be a positive integer. It must be one for
### univariate distributions and must be greater one for 
### multivariate distributions.
".init.valid.r.Model" <- function(obj)
{
    univ    <- .get.univ.Model()
    multiv  <- .get.multiv.Model()
    if (obj@r < 1) {
        warning(paste("Wrong specification of slot 'r' ",
                      "in 'model' object. Dimension of ",
                      "variables must be a positive integer.", 
                      sep =""))
    } else {
        if ((obj@dist %in% univ) && obj@r > 1) {
            warning(paste("Wrong specification of slot 'r' ",
                          "in 'model' object. Univariate ",
                          "distributions can only have one ",
                          "dimension.", sep = ""))
        } else if ((obj@dist %in% multiv) && obj@r < 2) {
            warning(paste("Wrong specification of slot 'r' ",
                          "in 'model' object. Multivariate ",
                          "distributions must have dimension ",
                          "greater one.", sep =""))
        }
    }
}

".valid.r.Model" <- function(obj)
{
    univ    <- .get.univ.Model()
    multiv  <- .get.multiv.Model() 
    if (obj@r < 1) {
        stop(paste("Wrong specification of slot 'r' ",
                   "in 'model' object. Dimension of ",
                   "variables must be positive.", 
                      sep =""))
    } else {
        if ((obj@dist %in% univ) && obj@r > 1) {
            stop(paste("Wrong specification of slot 'r' ",
                       "in 'model' object. Univariate ",
                       "distributions can only have one ",
                       "dimension.", sep = ""))
        } else if ((obj@dist %in% multiv) && obj@r < 2) {
            stop(paste("Wrong specification of slot 'r' ",
                       "in 'model' object. Multivariate ",
                       "distributions must have dimension ",
                       "greater one.", sep =""))
        }
    }
}

### Valid weight: The weights must be positive and add to one.
### Furthermore, @weight must be matrix of dimension 1 x K.
".init.valid.weight.Model" <- function(obj)
{
    if (!all(is.na(obj@weight))) {
        if (nrow(obj@weight) > 1) {
            warning(paste("Wrong dimension of slot 'weight' in ",
                          "'model' object. Dimension of slot ",
                          "'weight' must be 1 x K.", sep = ""))
        } else {
            if (ncol(obj@weight) != obj@K) {
                warning(paste("Wrong number of weights in slot 'weight' of ",
                              "'model' object. Number of weights does not ",
                              "match number of components in slot 'K'.", sep = ""))
            } else {
                if (is.integer(obj@weight)) {
                    stop(paste("Wrong specification of slot 'weight' of ",
                                  "'model' object. Weights must be of type ",
                                  "'numeric'.", sep = ""))
                }
                if (!is.numeric(obj@weight)) {
                     stop(paste("Wrong specification of slot 'weight' of ",
                                  "'model' object. Weights must be of type ",
                                  "'numeric'.", sep = ""))
                }
                if (any(obj@weight <= 0)) {
                    warning(paste("Weights in slot 'weight' of 'model' ",
                                  "object must be positive.", sep = ""))
                } else {
                    if (sum(obj@weight) != 1) {
                        warning(paste("Weights in slot 'weight' of 'model' ",
                                      "object must sum to one.", sep = ""))
                    }
                }
            }
        }
    }
}

".valid.weight.Model" <- function(obj)
{
    if (!all(is.na(obj@weight))) {
        if (nrow(obj@weight) > 1) {
            stop(paste("Wrong dimension of slot 'weight' in ",
                          "'model' object. Dimension of slot ",
                          "'weight' must be 1 x K.", sep = ""))
        } else {
            if (ncol(obj@weight) != obj@K) {
                stop(paste("Wrong number of weights in slot 'weight' of ",
                              "'model' object. Number of weights does not ",
                              "match number of components in slot 'K'.", sep = ""))
            } else {
                if (is.integer(obj@weight)) {
                    stop(paste("Wrong specification of slot 'weight' of ",
                                  "'model' object. Weights must be of type ",
                                  "'numeric'.", sep = ""))
                }
                if (!is.numeric(obj@weight)) {
                     stop(paste("Wrong specification of slot 'weight' of ",
                                  "'model' object. Weights must be of type ",
                                  "'numeric'.", sep = ""))
                }
                if (any(obj@weight <= 0)) {
                    stop(paste("Weights in slot 'weight' of 'model' ",
                                  "object must be positive.", sep = ""))
                } else {
                    if (sum(obj@weight) != 1) {
                        stop(paste("Weights in slot 'weight' of 'model' ",
                                      "object must sum to one.", sep = ""))
                    }
                }
            }
        }
    }
}

### Valid repetitions: Repetitions for Binomial models in @T must
### be a matrix of dimension N x 1.
### Repetitions must positive integers or NA.
".valid.T.Model" <- function(obj) 
{
    if (!all(is.na(obj@T))) {
        if (!is.integer(obj@T)) {
            stop(paste("Wrong type of slot 'T' in 'model' object ",
                       "Repetitions must be of type 'integer'.",
                       sep = ""))
        } 
        if (nrow(obj@T) > 1 && ncol(obj@T) > 1) {
            stop(paste("Wrong dimension of slot 'T' in 'model' ",
                       "object. Repetitions can only be ",
                       "one-dimensional", sep = ""))
        }
        if (any(obj@T < 1)) {
            stop(paste("Wrong specification of slot 'T' in 'model' ",
                       "object. Repetitions must be positive integers ",
                       "or NA.", sep = ""))
        }
    }
}

### Valid parameters: If @par is specified, parameters are checked
### for correcteness in regard to the specified model @dist.
### During initializing ('validObject()') the validity check 
### relies exclusively on errors, as it is assumed that the user
### knows what arguments have been provided to the constructor 
### 'model()'. The same is true for later usage of model objects.
### In case the user calls the setters of the model class it is 
### assumed, that the user intends to change the model object 
### slot by slot. In this case the validity check relies more on
### warnings.
".init.valid.par.Model" <- function(obj) 
{
    if (length(obj@par) > 0) {
        if (obj@dist %in% c("poisson", "cond.poisson")) {
            .init.valid.Poisson.Model(obj) 
        } else if (obj@dist == "binomial") {

        }
    }
}

".valid.par.Model" <- function(obj) 
{
    if (length(obj@par) > 0) {
        if (obj@dist %in% c("poisson", "cond.poisson")) {
            .valid.Poisson.Model(obj) 
        } else if (obj@dist == "binomial") {

        }
    }
}

### Valid Poisson: The Poisson parameters must be positive.
### Furthermore, the list element in @par must be named 
### 'lambda'. If model parameters are provided, the length
### of the object must be equal to @K.
### In case, that the user only modifies slots in the object
### it is relied on warnings instead of errors. Only in the
### case that parameters lambda in @par are specified and 
### are negative an error is thrown: Neither Poisson nor
### Exponential distributions allow negative- or zero-valued 
### parameters.
".init.valid.Poisson.Model" <- function(obj)
{
    if ("lambda" %in% names(obj@par)) {
        obj@par$lambda <- as.vector(obj@par$lambda)
        if (!is.numeric(obj@par$lambda) && !is.integer(obj@par$lambda)) {
            stop(paste("Wrong specification in slot 'par' of 'model' object. ",
                       "Parameters must be of type 'numeric' or 'integer'.",
                       sep = ""))
        }
        if (length(obj@par$lambda) != obj@K) {
            warning(paste("Wrong specification of slot 'par' of 'model' object. ", 
                          "Number of Poisson parameter in 'par$lambda' must ",
                          "match number of components in slot 'K'."), sep = "")
        } else {
            if (any(obj@par$lambda <= 0)) {
                stop(paste("Wrong specification of slot 'par' of 'model' ",
                              "object. Poisson parameters in 'par$lambda' ", 
                              "must be positive.", sep = ""))               
            }
        }
    } else {
        warning(paste("Wrong specification of slot 'par' in 'model' object. ",
                      "Poisson parameters must be named 'lambda'."), sep = "")
    }   
}

".valid.Poisson.Model" <- function(obj)
{
    if ("lambda" %in% names(obj@par)) {
        obj@par$lambda <- as.vector(obj@par$lambda)
        if (!is.numeric(obj@par$lambda) && !is.integer(obj@par$lambda)) {
            stop(paste("Wrong specification in slot 'par' of 'model' object. ",
                       "Parameters must be of type 'numeric' or 'integer'.",
                       sep = ""))
        }
        if (length(obj@par$lambda) != obj@K) {
            stop(paste("Wrong specification of slot 'par' of 'model' object. ", 
                       "Number of Poisson parameter in 'par$lambda' must ",
                       "match number of components in slot 'K'."), sep = "")
        } else {
            if (any(obj@par$lambda <= 0)) {
                stop(paste("Wrong specification of slot 'par' of 'model' ",
                           "object. Poisson parameters in 'par$lambda' ", 
                           "must be positive.", sep = ""))               
            }
        }
    } else {
        stop(paste("Wrong specification of slot 'par' in 'model' object. ",
                   "Poisson parameters must be named 'lambda'."), sep = "")
    }   
}

### valid.Binomial
".valid.Binomial.Model" <- function(model.obj) 
{
    if (dim(model.obj@T)[1] > 1 && dim(model.obj@T)[2] > 1) {
        stop(paste("Dimensions of repetitions 'T' for binomial mixture", 
             "model do not match conditions. Only one-dimensional",
             "repetitions can be used in a binomial mixture model."), sep ="")
    }
    
}

### Additional functions
".get.univ.Model" <- function()
{
    univ <- c("poisson", "cond.poisson",
              "binomial", "exponential",
              "normal", "student")
    return(univ)
}

".get.multiv.Model" <- function()
{
    multiv <- c("normult", "studmult")
    return(multiv)
}
