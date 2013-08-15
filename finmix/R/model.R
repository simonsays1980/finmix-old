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
                       dist <- object@dist
                       .valid.Dist(object)
                       .valid.Weight(object)
                       if (object@K <= 0) {					
                           stop("Number of components 'K' must be a 
                                positive integer.")
                       }
                       if (dist == "poisson") {
                           object <- .valid.Poisson.Model(object)
                       } else if (dist == "binomial") {
                           .valid.Binomial.Model(object)
                       } else {

                       }
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
"model" <- function(dist = "poisson", r = as.integer(1), 
                    K = as.integer(1), weight = matrix(), 
                    par = list(), indicmod = "multinomial", 
           			indicfix = FALSE, T = matrix()) 
{
    if (K > 1 && all(is.na(weight))) {
	    weight <- matrix(1/K, nrow = 1, ncol = K)
	} 	
    object <- .model(dist = dist, r = as.integer(r), 
                     K = as.integer(K), weight = weight, par = par, 
                     indicmod = indicmod, indicfix = indicfix, T = T)
    return(object)
}

setMethod("plot", "model", 
          function(x, y, ..., dev = TRUE) 
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
              }	else if (dist == "poisson" || dist == "cond.poisson") {
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

## Marginal Mixture ##
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))
setMethod("mixturemar", "model", 
          function(object, J) 
          {
              if (object@dist == "normult") {	
                  dist <- ifelse(length(J) == 1, "normal", "normult")
                  r             <- length(J)
                  K             <- object@K
                  weight        <- object@weight
                  mu            <- object@par$mu[J, ]
                  sigma         <- object@par$sigma[J, J, ]
                  par           <- list(mu = mu, sigma = sigma)
                  indicmod      <- "multinomial"
                  indicfix      <- TRUE
                  margin.model  <- .model(dist = dist, r = r, K = K, 
                                          weight = weight, par = par, 
                                          indicmod = indicmod, 
                                          indicfix = indicfix)
                  validObject(margin.model)
                  return(margin.model)
              } else if (object@dist == "studmult") {
                  dist <- ifelse(length(J) == 1, "student", "studmult")
                  r             <- length(J)
                  K             <- object@K
                  weight        <- object@weight
                  mu            <- object@par$mu[J, ]
                  sigma         <- object@par$sigma[J, J, ] 
                  df            <- object@par$df
                  par           <- list(mu = mu, sigma = sigma, df = df)
                  indicmod      <- "multinomial"
                  indicfix      <- TRUE
                  margin.model  <- .model(dist = dist, r = r, K = K, 
                                          weight = weight, par = par, 
                                          indicmod = indicmod,
                                          indicfix = indicfix)
                  validObject(margin.model)
                  return(margin.model)
              } else {
                  stop("A marginal distribution can only be obtained from 
                       multivariate distributions.")
              }
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
              if (!all(is.na(object@par))) {
                 cat("     weight      :", 
                      paste(dim(object@weight), collapse = "x"), 
                      "\n")
                 cat("     par         : List of",
                      length(object@par), "\n")
              }
              cat("     indicmod    :", object@indicmod, "\n")
              cat("     indicfix    :", object@indicfix, "\n")
              if (object@dist == "binomial") {
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
                     validObject(object)
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
                     object@K       <- as.integer(value)
                     object@weight  <- matrix(1/value, 
                                              nrow = 1, 
                                              ncol = object@K)
                     validObject(object)
                     return(object)
                 }
)

setReplaceMethod("setWeight", "model", 
                 function(object, value) 
                 {
                     object@weight <- value
                     validObject(object)
                     return(object)
                 }
)

setReplaceMethod("setPar", "model", 
                 function(object, value) 
                 {
                     object@par <- value
                     validObject(object)
                     return(object)
                 }
)

setReplaceMethod("setIndicmod", "model", 
                 function(object, value) 
                 {
                     object@indicmod <- value
                     validObject(object)
                     return(object)
                 }
)

setReplaceMethod("setIndicfix", "model", 
                 function(object, value) 
                 {
                     object@indicfix <- value
                     validObject(object)
                     return(object)
                 }
)

setReplaceMethod("setT", "model",
                 function(object, value) 
                 {
                     object@T <- as.matrix(value)
                     validObject(object)
                     return(object)
                 }
)

### Private functions
### These functions are not exported
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

### Validity
### Valid.dist: @dist must be a distribution which
### is implemented
".valid.Dist" <- function(model.obj)
{
    dists            <- c("normal", "normult", "exponential", 
                          "student", "studmult", "poisson", 
                          "cond.poisson", "binomial")
    indicmod.dists   <- c("multinomial")
    if (length(model.obj@dist) > 0) {
        if (!(model.obj@dist %in% dists)) {
            stop("Unknown distribution in slot 'dist'.")
        } else {
            if (!(model.obj@indicmod) %in% indicmod.dists) {
            stop(paste("Unknown indicator distribution in slot", 
                 "'indicmod'"), sep = "")
            }
        }
    }
}

### valid.Weight: The weights must be positive and add to one.
".valid.Weight" <- function(model.obj)
{
    if (!all(is.na(model.obj@weight))) {
        if (ncol(model.obj@weight) != model.obj@K) {
            warning(paste("Number of weights in slot 'weight' do not", 
                    "match number of components in slot 'K'."), sep = "")
        } else {
            if (any(model.obj@weight <= 0)) {
                stop("Weights in slot 'weight' must be positive.")
            } else {
                if (sum(model.obj@weight) != 1) {
                    stop("Weights in slot 'weight' must sum to one.")
                }
            }
        }
    }
}
### valid.Poisson: The Poisson parameters must be positive.
### Furthermore, the list element in @par must be named 
### 'lambda'. If model parameters are provided, the length
### of the object must be equal to @K.
".valid.Poisson.Model" <- function(model.obj)
{
    if (length(model.obj@par) > 0) {
        if ("lambda" %in% names(model.obj@par)) {
            model.obj@par$lambda <- as.vector(model.obj@par$lambda)
            if (length(model.obj@par$lambda) != model.obj@K) {
                warning(paste("Number of Poisson parameters in slot 'par'", 
                        "do not match number of components in slot", 
                        "'K'."), sep = "")
            } else {
                if (any(model.obj@par$lambda <= 0)) {
                    stop(paste("Poisson parameters 'lambda' in slot", 
                            "'par' must be positive."), sep = "")
                }
            }
        } else {
            warning(paste("Poisson parameters in slot 'par' must be", 
                    "named 'lambda'."), sep = "")
        }
    }
    return(model.obj)
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
