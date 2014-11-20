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

.mcmcoutputpermhier <- setClass( "mcmcoutputpermhier",
                                 contains = c( "mcmcpermind", 
                                               "mcmcoutputhier" ),
                                 validity = function( object ) 
                                 {
                                     ## else: OK
                                     TRUE
                                 }
)

setMethod("initialize", "mcmcoutputpermhier",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), relabel = character(), 
                   weightperm = array(), logperm = list(), 
                   entropyperm = array(), STperm = array(), 
                   Sperm = array(), NKperm = array()) 
          {
              .Object@M             <- mcmcoutput@M
              .Object@burnin        <- mcmcoutput@burnin
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@hyper         <- mcmcoutput@hyper
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@relabel       <- relabel
              .Object@weightperm    <- weightperm
              .Object@logperm       <- logperm
              .Object@entropyperm   <- entropyperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermhier", 
          function(object)
          {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     relabel     :", object@relabel, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              if (!all(is.na(object@S))) {
                  cat("     S           :", 
                      paste(dim(object@S), collapse = "x"), "\n")
              }
              cat("     NK          :",
                  paste(dim(object@NK), collapse = "x"), "\n")
              cat("     clust       :",
                  paste(dim(object@clust), collapse = "x"), "\n")
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     weightperm  :",
                  paste(dim(object@weightperm), collapse = "x"), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     entropyperm :",
                  paste(dim(object@entropyperm), collapse = "x"), "\n")
              cat("     STperm      :",
                  paste(dim(object@STperm), collapse = "x"), "\n")
              if (!all(is.na(object@Sperm))) {
                  cat("     Sperm       :",
                      paste(dim(object@Sperm), collapse = "x"), "\n")
              }
              cat("     NKperm      :", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputpermhier",
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
          {
              dist <- x@model@dist
              if ( lik %in% c( 0, 1 ) ) {
                  if ( dist == "poisson" ) {
                      .permtraces.Poisson.Base.Hier( x, dev )
                  } else if (dist == "binomial") {
                      .permtraces.Binomial.Base( x, dev )
                  } else if ( dist == "exponential" ) {
                      .permtraces.Exponential.Base( x, dev )
                  } else if ( dist == "normal" ) {
                      .permtraces.Normal.Hier( x, dev )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "student" ) {
                      .permtraces.Student.Hier( x, dev )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "normult" ) {
                      .permtraces.Normult.Hier( x, dev, col )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "studmult" ) {
                      .permtraces.Studmult.Hier( x, dev, col )
                      .permtraces.Weights.Base( x, dev, col )
                  }
              }	
              if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .permtraces.Log.Base( x, dev, col )
              }
          }
)

setMethod( "plotHist", signature( x = "mcmcoutputpermhier", 
                                  dev = "ANY" ), 
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if( dist == "poisson" ) {
                   .permhist.Poisson.Base.Hier( x, dev )
               } else if ( dist == "binomial" ) {
                   .permhist.Binomial.Base( x, dev )
               } else if ( dist == "exponential" ) {
                   .permhist.Exponential.Base( x, dev )
               } else if ( dist == "normal" ) {
                   .permhist.Normal.Base.Hier( x, dev )
               } else if ( dist == "student" ) {
                   .permhist.Student.Base.Hier( x, dev )
               } else if ( dist == "normult" ) {
                   .permhist.Normult.Base.Hier( x, dev )
               } else if ( dist == "studmult"  ) {
                   .permhist.Studmult.Base.Hier( x, dev )
               }
           }
)

setMethod( "plotDens", signature( x = "mcmcoutputpermhier", 
                                  dev = "ANY" ), 
           function(x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permdens.Poisson.Base.Hier( x, dev )
               } else if ( dist == "binomial" ) {
                   .permdens.Binomial.Base( x, dev )
               } else if ( dist == "exponential" ) {
                   .permdens.Exponential.Base( x, dev )
               } else if ( dist == "normal" ) {
                   .permdens.Normal.Base.Hier( x, dev )
               } else if ( dist == "student" ) {
                   .permdens.Student.Base.Hier( x, dev )
               } else if ( dist == "normult" ) {
                   .permdens.Normult.Base.Hier( x, dev )
               } else if ( dist == "studmult"  ) {
                   .permdens.Studmult.Base.Hier( x, dev )
               }
 
            }
)

setMethod( "plotPointProc", signature( x      = "mcmcoutputpermhier",
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

setMethod( "plotSampRep", signature( x    = "mcmcoutputpermhier",
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

setMethod( "plotPostDens", signature( x   = "mcmcoutputpermhier",
                                      dev = "ANY" ),
           function(x, dev = TRUE, ... ) 
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

### Plot
### Traces
### Traces Poisson: Plots the traces for all Poisson 
### parameters, the weights and the hpyer-parameter 'b'.
".permtraces.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K
		trace.n <- K * 2
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots (permuted)")
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
    weight <- x@weightperm
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
}

### Histograms
### Histograms Poisson: plots histograms for all Poisson parameters,
### the weights and the hyper-parameter 'b'.
".permhist.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    lambda <- x@parperm$lambda
    weight <- x@weightperm
    b <- x@hyper$b
    vars        <- cbind(lambda, weight[, seq(1:(K - 1))], b)
    lab.names   <- vector("list", 2 * K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in (K + 1):(2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }
    lab.names[[2 * K]] <- "b"
    .symmetric.Hist(vars, lab.names)
}

".permhist.Normal.Base.Hier" <- function( x, dev )
{
    .permhist.Normal( x, dev )
    C   <- x@hyper$C
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histogram Hyperparameter" )
    }    
    .symmetric.Hist( C, "C" )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weights.lab.names   <- vector( "list", K ) 
        for ( k in 1:K ) {
            weights.lab.names[[k]]  <- bquote( eta[.( k )] )
        }
        if ( K > 1 ) {
            if ( .check.grDevice() && dev ) {
                dev.new( title = "Histograms Weights (permuted)" )
            }
            .symmetric.Hist( weight, weights.lab.names )
        }
    }
}

".permhist.Student.Base.Hier" <- function( x, dev )
{
    .permhist.Student( x, dev )
    C   <- x@hyper$C
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histogram Hyperparameter" )
    }    
    .symmetric.Hist( C, "C" )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weight.lab.names    <- vector( "list", K )
        for ( k in 1:K ) {
            weight.lab.names[[k]]   <- bquote( eta[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Weights (permuted)" )
        }
        .symmetric.Hist( weight, weight.lab.names )
    }
}

".permhist.Normult.Base.Hier"  <- function( x, dev ) 
{
    .permhist.Normult( x, dev )   
    Clogdet <- sapply( seq( 1:x@M ), function( i ) log( det( qinmatr( x@hyper$C[i,] ) ) ) )
    Ctr     <- sapply( seq( 1:x@M ), function( i ) sum( diag( qinmatr( x@hyperC[i,]) ) ) )
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histograms Hyperparameters" )        
    }
    lab.C.names <- vector( "list", 2 )
    lab.C.names[[1]]    <- "tr(C)"
    lab.C.names[[2]]    <- "log(det(C))"
    vars                <- cbind( Ctr, Clogdet )
    .symmetric.Hist( vars, lab.C.names )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weight.lab.names    <- vector( "list", K )        
        for ( k in 1:K ) {
            weight.lab.names[[k]]   <- bquote( eta[.( k ) ] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Weights (permuted)" )
        }
        .symmetric.Hist( weight, weight.lab.names )
    }
}

".permhist.Studmult.Base.Hier"  <- function( x, dev ) 
{
    .permhist.Studmult( x, dev )
    Clogdet <- sapply( seq( 1:x@M ), function( i ) log( det( qinmatr( x@hyper$C[i,] ) ) ) )
    Ctr     <- sapply( seq( 1:x@M ), function( i ) sum( diag( qinmatr( x@hyperC[i,]) ) ) )
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histograms Hyperparameters" )        
    }
    lab.C.names <- vector( "list", 2 )
    lab.C.names[[1]]    <- "tr(C)"
    lab.C.names[[2]]    <- "log(det(C))"
    vars                <- cbind( Ctr, Clogdet )
    .symmetric.Hist( vars, lab.C.names )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weight.lab.names    <- vector( "list", K )
        for ( k in 1:K ) {
            weight.lab.names[[k]]   <- bquote( eta[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Weights (permuted)" )
        }
       .symmetric.Hist( weight, weight.lab.names )
    }
}

### Densities
### Densities Poisson: plots Kernel densities for all Poisson 
### parameters, the weights and the hyper-parameter 'b'.
".permdens.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    lambda <- x@parperm$lambda
    weight <- x@weightperm
    b <- x@hyper$b
    vars        <- cbind(lambda, weight[, seq(1:(K - 1))], b)
    lab.names   <- vector("list", 2 * K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in (K + 1):(2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }
    lab.names[[2 * K]] <- "b"
    .symmetric.Dens(vars, lab.names)
}

".permdens.Normal.Base.Hier" <- function( x, dev )
{
    .permdens.Normal( x, dev )
    C   <- x@hyper$C
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Density Hyperparameter" )
    }    
    .symmetric.Dens( C, "C" )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weights.lab.names   <- vector( "list", K ) 
        for ( k in 1:K ) {
            weights.lab.names[[k]]  <- bquote( eta[.( k )] )
        }
        if ( K > 1 ) {
            if ( .check.grDevice() && dev ) {
                dev.new( title = "Densities Weights (permuted)" )
            }
            .symmetric.Dens( weight, weights.lab.names )
        }
    }
}

".permdens.Student.Base.Hier" <- function( x, dev )
{
    .permdens.Student( x, dev )
    C   <- x@hyper$C
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histogram Hyperparameter" )
    }    
    .symmetric.dens( C, "C" )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weight.lab.names    <- vector( "list", K )
        for ( k in 1:K ) {
            weight.lab.names[[k]]   <- bquote( eta[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Weights (permuted)" )
        }
        .symmetric.Dens( weight, weight.lab.names )
    }
}

".permdens.Normult.Base.Hier"  <- function( x, dev ) 
{
    .permdens.Normult( x, dev )   
    Clogdet <- sapply( seq( 1:x@M ), function( i ) log( det( qinmatr( x@hyper$C[i,] ) ) ) )
    Ctr     <- sapply( seq( 1:x@M ), function( i ) sum( diag( qinmatr( x@hyperC[i,]) ) ) )
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Densities Hyperparameters" )        
    }
    lab.C.names <- vector( "list", 2 )
    lab.C.names[[1]]    <- "tr(C)"
    lab.C.names[[2]]    <- "log(det(C))"
    vars                <- cbind( Ctr, Clogdet )
    .symmetric.Dens( vars, lab.C.names )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weight.lab.names    <- vector( "list", K )        
        for ( k in 1:K ) {
            weight.lab.names[[k]]   <- bquote( eta[.( k ) ] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Weights (permuted)" )
        }
        .symmetric.Dens( weight, weight.lab.names )
    }
}

".permdens.Studmult.Base.Hier"  <- function( x, dev ) 
{
    .permdens.Studmult( x, dev )
    Clogdet <- sapply( seq( 1:x@M ), function( i ) log( det( qinmatr( x@hyper$C[i,] ) ) ) )
    Ctr     <- sapply( seq( 1:x@M ), function( i ) sum( diag( qinmatr( x@hyperC[i,]) ) ) )
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Densities Hyperparameters" )        
    }
    lab.C.names <- vector( "list", 2 )
    lab.C.names[[1]]    <- "tr(C)"
    lab.C.names[[2]]    <- "log(det(C))"
    vars                <- cbind( Ctr, Clogdet )
    .symmetric.Dens( vars, lab.C.names )
    if ( K > 1 ) {
        weight              <- x@weightperm
        weight.lab.names    <- vector( "list", K )
        for ( k in 1:K ) {
            weight.lab.names[[k]]   <- bquote( eta[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Weights (permuted)" )
        }
       .symmetric.Dens( weight, weight.lab.names )
    }
}

