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

.mcmcoutputpermfixhier <- setClass( "mcmcoutputpermfixhier",
                                    contains = c( "mcmcpermfix", "mcmcoutputfixhier" ),
                                    validity = function( object ) 
                                    {
                                        ## else: OK
                                        TRUE
                                    }
)

setMethod("initialize", "mcmcoutputpermfixhier",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), logperm = list()) 
          {
              .Object@M         <- mcmcoutput@M
              .Object@burnin    <- mcmcoutput@burnin
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@hyper     <- mcmcoutput@hyper
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object@logperm   <- logperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixhier",
          function(object) 
          {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputpermfixhier", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
          {
              dist <- x@model@dist
              if ( lik %in% c( 0, 1 ) ) {                  
                  if ( dist == "poisson" ) {
                      .permtraces.Poisson.Hier( x, dev )
                  } else if ( dist == "binomial" ) {
                      .permtraces.Binomial( x, dev )
                  } else if ( dist == "exponential" ) {
                      callNextMethod( x, dev, lik, col, ... )
                  } else if ( dist == "normal" ) {
                      .permtraces.Normal.Hier( x, dev )
                  } else if ( dist == "student" ) {
                      .permtraces.Student.Hier( x, dev )
                  } else if ( dist == "normult" ) {
                      .permtraces.Normult.Hier( x, dev, col )
                  } else if ( dist == "studmult" ) {
                      .permtraces.Studmult.Hier( x, dev, col )
                  }
              }
              if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .permtraces.Log( x, dev, col )
              }
          }
)

setMethod( "plotHist", signature( x   = "mcmcoutputpermfixhier", 
                                  dev = "ANY" ), 
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permhist.Poisson.Hier( x, dev )
               } else if ( dist == "binomial" ) {
                   .permhist.Binomial( x, dev )
               } else if ( dist == "exponential" ) {
                   .permhist.Exponential( x, dev )
               } else if ( dist == "normal" ) {
                   .permhist.Normal.Hier( x, dev )
               } else if ( dist == "student" ) {
                   .permhist.Student.Hier( x, dev )
               } else if ( dist == "normult" ) {
                   .permhist.Normult.Hier( x, dev )
               } else if ( dist == "studmult" ) {
                   .permhist.Studmult.Hier( x, dev )
               }
           }
)

setMethod( "plotDens", signature( x   = "mcmcoutputpermfixhier", 
                                  dev = "ANY" ), 
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permdens.Poisson.Hier( x, dev )
               } else if ( dist == "binomial" ) {
                   .permdens.Binomial( x, dev )
               } else if ( dist == "exponential" ) {
                   .permdens.Exponential( x, dev )
               } else if ( dist == "normal" ) {
                   .permdens.Normal.Hier( x, dev )
               } else if ( dist == "student" ) {
                   .permdens.Student.Hier( x, dev )
               } else if ( dist == "normult" ) {
                   .permdens.Normult.Hier( x, dev )
               } else if ( dist == "studmult" ) {
                   .permdens.Studmult.Hier( x, dev )
               }

           }
)

setMethod( "plotPointProc", signature( x      = "mcmcoutputpermfixhier",
                                      dev    = "ANY" ),
           function( x, dev = TRUE, ... )
           {
               dist <- x@model@dist
               if ( dist %in% c( "poisson", "exponential" ) ) {
                   .permpointproc.Poisson( x, dev )
               } else if ( dist == "binomial" ) {
                   .permpointproc.Binomial( x, dev )
               } else if ( dist == "exponential" ) {
                   .permpointproc.Exponential( x, dev )
               } else if ( dist %in% c( "normal", "student" ) ) {
                   .permpointproc.Normal( x, dev )
               } else if ( dist %in% c( "normult", "studmult" ) ) {
                   .permpointproc.Normult( x, dev )
               }
           }
)

setMethod( "plotSampRep", signature( x    = "mcmcoutputpermfixhier",
                                     dev  = "ANY" ),
           function( x, dev, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permsamprep.Poisson( x, dev )
               } else if ( dist == "binomial" ) {
                   .permsamprep.Binomial( x, dev )                                      
               } else if ( dist == "exponential" ) {
                   .permsamprep.Exponential( x, dev )
               } else if ( dist == "normal" ) {
                   .permsamprep.Normal( x, dev )
               } else if ( dist == "student" ) {
                   .permsamprep( x, dev )
               } else if ( dist == "normult" ) {
                   .permsamprep.Normal( x, dev)
               } else if ( dist == "studmult" ) {
                   .permsamprep.Studmult( x, dev )
               }
           }
)

setMethod( "plotPostDens", signature( x   = "mcmcoutputpermfixhier",
                                      dev = "ANY" ),
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist %in% c( "poisson", "exponential" ) ) {
                   .permpostdens.Poisson( x, dev )
               } else if ( dist == "binomial" ) {
                   .permpostdens.Binomial( x, dev )
               } else if ( dist == "normal" ) {
                   .permpostdens.Normal( x, dev )
               } else if ( dist == "student" ) {
                   .permpostdens.Student( x, dev )
               } else if ( dist == "normult" ){
                   .permpostdens.Normult( x, dev )
               } else if ( dist == "studmult" ) {
                   .permpostdens.Studmult( x, dev )
               }
           }
)

### Private functions.
### These functions are not exported.

### Traces 
### Traces Poisson: Plots the traces of Poisson parameters 
### and the hyper-parameter 'b'.
".permtraces.Poisson.Hier" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K + 1
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    lambda <- x@parperm$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]),
              cex = 0.6, line = 3)
    }
    b <- x@hyper$b
    plot(b, type = "l", axes = F, 
         col = "gray68", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

".permtraces.Normal.Hier"   <- function( x, dev ) 
{
    K       <- x@model@K
    trace.n <- 2 * K + 1
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" )        
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    for ( k in 1:K ) {
        plot( mu[, k], type = "l", axes = F,
             col = "gray20", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu[k = .( k )] ),
               cex = .6, line = 3 )       
    }
    for ( k in 1:K ) {
        plot( sigma[, k], type = "l", axes = F,
              col = "gray30", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )              
        mtext( side = 2, las = 2, bquote( sigma[k = .( k )]),
               cex = .6, line = 3 )        
    }
    C   <- x@hyper$C
    plot( c, type = "l", axes = F, 
          col = "gray68", xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = 0.7 )
    mtext( side = 2, las = 2, "C", cex = .6, line = 3 )
    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )   
}

### --------------------------------------------------------------------
### .permtraces.Student.Hier
### @description    Plots traces for parameters of a univariate Student 
###                 mixture.
### @par    x       an object of class mcmcoutputfix
###         dev     an object of class 'logical' 
### @detail         Plots the traces for each component parameter of an
###                 Student mixture. If 'dev' is set to FALSE 
###                 (TRUE is default) no device is created, instead 
###                 the graphic can be stored to a file.
### @see            ?mcmcoutput, ?plotTraces
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------
".permtraces.Student.Hier" <- function( x, dev ) 
{
    K       <- x@model@K
    trace.n <- 3 * K + 1
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" )        
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    df      <- x@parperm$df
    for ( k in 1:K ) {
        plot( mu[, k], type = "l", axes = F,
             col = "gray20", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu[k = .( k )] ),
               cex = .6, line = 3 )       
    }
    for ( k in 1:K ) {
        plot( sigma[, k], type = "l", axes = F,
              col = "gray30", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )              
        mtext( side = 2, las = 2, bquote( sigma[k = .( k )]),
               cex = .6, line = 3 )        
    }
    for ( k in 1:K ) {
        plot( df[, k], type = "l", axes = F,
              col = "gray40", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( nu[k = .( k )]),
               cex = .6, line = 3 )
    }
    C <- x@hyper$C
    plot( C, type = "l", axes = F,
          col = "gray68", xlab = "", ylab = "" ) 
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, "C", cex = .6,
           line = 3 )
    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
}

"permtraces.Normult.Hier" <- function( x, dev, col ) 
{
    .permtraces.Normult( x, dev, col )
    r   <- x@model@r
    K   <- x@model@K
    C   <- x@hyper$C 
    C.trace     <- sapply( seq( 1, x@M ),
                           function( i ) sum( diag( qinmatr( C[i,] ) ) ) )
    C.logdet    <- sapply( seq( 1, x@M ),
                           function( i ) log( det( qinmatr( C[i,] ) ) ) )
    # C traces
    mmax    <- max( C.trace )
    mmin    <- min( C.trace )
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots Hyperparameters" )
    }
    par( mfrow = c( 2, 1 ), mar = c( 1, 2, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    if ( col ) {
        cscale  <- rainbow( K, start = 0.5, end = 0 )        
    } else {
        cscale  <- gray.colors( K, start = 0.5, end = 0.15 ) 
    }
    plot( C.trace, type = "l", axes = F,
          col = cscale[K], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( tr(C) ),
           cex = .6, line = 3 )
    plot( C.logdet, type = "l", axes = F,
          col = cscale[K], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = .7 )
    name    <- vector( "character", K )
    mtext( side = 2, las = 2, bquote( log(det(C))),
           cex = .6, line = 3 )
    axis( 1 )
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
}

".permtraces.Studmult.Hier" <- function( x, dev, col ) 
{
    .permtraces.Studmult( x, dev, col )
    r   <- x@model@r
    K   <- x@model@K
    C   <- x@hyper$C 
    C.trace     <- sapply( seq( 1, x@M ),
                           function( i ) sum( diag( qinmatr( C[i,] ) ) ) )
    C.logdet    <- sapply( seq( 1, x@M ),
                           function( i ) log( det( qinmatr( C[i,] ) ) ) )
    
    # C traces
    mmax    <- max( C.trace )
    mmin    <- min( C.trace )
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots Hyperparameters" )
    }
    par( mfrow = c( 2, 1 ), mar = c( 1, 2, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    if ( col ) {
        cscale  <- rainbow( K, start = 0.5, end = 0 )        
    } else {
        cscale  <- gray.colors( K, start = 0.5, end = 0.15 ) 
    }
    plot( C.trace, type = "l", axes = F,
          col = cscale[K], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( tr(C) ),
           cex = .6, line = 3 )
    plot( C.logdet, type = "l", axes = F,
          col = cscale[K], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( log(det(C))),
           cex = .6, line = 3 )
    axis( 1 )
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
}

### Histograms
### Histograms Poisson: Plots histograms for all Poisson 
### parameters and the hyper-parameter 'b'.
".permhist.Poisson.Hier." <- function(x, dev)
{
    K <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
	}
    lambda  <- x@parperm$lambda
    b       <- x@hyper$b
    vars    <- cbind(lambda, b)
    lab.names <- vector("list", K + 1)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    lab.names[[K + 1]] <- "b"
    .symmetric.Hist(vars, lab.names)
}

".permhist.Normal.Hier" <- function( x, dev )
{
    K       <- x@model@K 
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    C       <- x@hyper$C
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Mu (permuted)" )
        }
        .symmetric.Hist( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram (permuted)" )
        }
        .symmetric.Hist( sigma, list( bquote( sigma ) ) ) 
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Mu (permuted)" )
        }
        .symmetric.Hist( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Sigma (permuted)" )
        }
        .symmetric.Hist( sigma, sigma.lab.names )
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histogram Hyperparameter C" )
    }
    .symmetric.Hist( C, "C" )
}

".permhist.Student.Hier" <- function( x, dev )
{
    K       <- x@model@K 
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    degf    <- x@parperm$df
    C       <- x@hyper$C
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Mu (permuted)" )
        }
        .symmetric.Hist( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Sigma (permuted)" )
        }
        .symmetric.Hist( sigma, list( bquote( sigma ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Degrees of Freedom (permuted)" )
        }
        .symmetric.Hist( degf, list( bquote( nu ) ) )
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            degf.lab.names[[k]]     <- bquote( nu[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Mu (permuted)" )
        }
        .symmetric.Hist( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Sigma (permuted)" )
        }
        .symmetric.Hist( sigma, sigma.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Degrees of Freedom (permuted)" )
        }
        .symmetric.Hist( degf, degf.lab.names )
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histogram Hyperparameter C" )
    }
    .symmetric.Hist( C, "C" )
}

".permhist.Normult.Hier"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    logdetC <- sapply( seq( 1, x@M ), function( i ) log( det( qinmatr( x@hyper$C[i, ] ) ) ) ) 
    trC     <- sapply( seq( 1, x@M ), function( i ) sum( diag( qinmatr( x@hyper$C[i, ] ) ) ) )
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], list( bquote( sigma ) ) )           
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], sigma.lab.names )
        }        
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histograms Hyperparameter C" ) 
    } 
    C.lab.names         <- vector( "list", 2 )
    C.lab.names[[1]]    <- "log(det(C))"
    C.lab.names[[2]]    <- "tr(C)"    
    .symmetric.Hist( cbind( logdetC, trC ), C.lab.names )
}

".permhist.Studmult.Hier"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    degf    <- x@parperm$df
    logdetC <- sapply( seq( 1, x@M ), function( i ) log( det( qinmatr( x@hyper$C[i, ] ) ) ) )
    trC     <- sapply( seq( 1, x@M ), function( i ) sum( diag( qinmatr( x@hyper$C[i, ] ) ) ) )
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], list( bquote( sigma ) ) )
         
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], sigma.lab.names )
        }
    }
    if ( K == 1 ) { 
        if (.check.grDevice() & dev ) {
            dev.new( title = paste( "Histograms Feature ", rr, 
                                    " Mu (permuted)", sep = "" ) )
        }
        .symmetric.Hist( degf[, rr,], list( bquote( nu ) ) ) 
    } else {
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            degf.lab.names[[k]] <- bquote( nu[.( k )] )
        }
        if (.check.grDevice() & dev ) {
            dev.new( title = paste( "Histograms Feature ", rr, 
                                    " Sigma (permuted)", sep = "" ) )
        }
        .symmetric.Hist( degf[, rr,], degf.lab.names ) 
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histograms Hyperparameter C" ) 
    } 
    C.lab.names         <- vector( "list", 2 )
    C.lab.names[[1]]    <- "log(det(C))"
    C.lab.names[[2]]    <- "tr(C)"    
    .symmetric.Hist( cbind( logdetC, trC ), C.lab.names )

}

### Densities
### Densities Poisson: Plots densities for all Poisson 
### parameters and the hyper-parameter 'b'.
".permdens.Poisson.Hier." <- function(x, dev)
{
    K <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
	}
    lambda  <- x@parperm$lambda
    b       <- x@hyper$b
    vars    <- cbind(lambda, b)
    lab.names <- vector("list", K + 1)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    lab.names[[K + 1]] <- "b"
    .symmetric.Dens(vars, lab.names)
}

".permdens.Normal.Hier" <- function( x, dev )
{
    K       <- x@model@K 
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    C       <- x@hyper$C
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Mu (permuted)" )
        }
        .symmetric.Dens( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Sigma (permuted)" )
        }
        .symmetric.Dens( sigma, list( bquote( sigma ) ) ) 
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Mu (permuted)" )
        }
        .symmetric.Dens( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Sigma (permuted)" )
        }
        .symmetric.Dens( sigma, sigma.lab.names )
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histogram Hyperparameter C" )
    }
    .symmetric.Dens( C, "C" )
}

".permdens.Student.Hier" <- function( x, dev )
{
    K <- x@model@K 
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    degf    <- x@parperm$df
    C       <- x@hyper$C
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Mu (permuted)" )
        }
        .symmetric.Dens( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Sigma (permuted)" )
        }
        .symmetric.Dens( sigma, list( bquote( sigma ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Degrees of Freedom (permuted)" )
        }
        .symmetric.Dens( degf, list( bquote( nu ) ) )
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            degf.lab.names[[k]]     <- bquote( nu[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Mu (permuted)" )
        }
        .symmetric.Dens( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Sigma (permuted)" )
        }
        .symmetric.Dens( sigma, sigma.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Degrees of Freedom (permuted)" )
        }
        .symmetric.Dens( degf, degf.lab.names )
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Density Hyperparameter C" )
    }
    .symmetric.Dens( C, "C" )
}

".permdens.Normult.Hier"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    logdetC <- sapply( seq( 1, x@M ), function( i ) log( det( qinmatr( x@hyper$C[i, ] ) ) ) ) 
    trC     <- sapply( seq( 1, x@M ), function( i ) sum( diag( qinmatr( x@hyper$C[i, ] ) ) ) )
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], list( bquote( sigma ) ) )           
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], sigma.lab.names )
        }        
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Densities Hyperparameter C" ) 
    } 
    C.lab.names         <- vector( "list", 2 )
    C.lab.names[[1]]    <- "log(det(C))"
    C.lab.names[[2]]    <- "tr(C)"    
    .symmetric.Dens( cbind( logdetC, trC ), C.lab.names )
}

".permdens.Studmult.Hier"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
    degf    <- x@parperm$df
    logdetC <- sapply( seq( 1, x@M ), function( i ) log( det( qinmatr( x@hyper$C[i, ] ) ) ) )
    trC     <- sapply( seq( 1, x@M ), function( i ) sum( diag( qinmatr( x@hyper$C[i, ] ) ) ) )
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], list( bquote( sigma ) ) )
         
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu (permuted)", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma (permuted)", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], sigma.lab.names )
        }
    }
    if ( K == 1 ) { 
        if (.check.grDevice() & dev ) {
            dev.new( title = "Density Degrees of Freedom (permuted)" )                                      
        }
        .symmetric.Dens( degf[, rr,], list( bquote( nu ) ) ) 
    } else {
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            degf.lab.names[[k]] <- bquote( nu[.( k )] )
        }
        if (.check.grDevice() & dev ) {
            dev.new( title = "Densities Degrees of Freedom (permuted)" )                                    
        }
        .symmetric.Dens( degf[, rr,], degf.lab.names ) 
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Densities Hyperparameter C" ) 
    } 
    C.lab.names         <- vector( "list", 2 )
    C.lab.names[[1]]    <- "log(det(C))"
    C.lab.names[[2]]    <- "tr(C)"    
    .symmetric.Dens( cbind( logdetC, trC ), C.lab.names )

}

