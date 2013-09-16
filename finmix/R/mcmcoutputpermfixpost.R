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

.mcmcoutputpermfixpost <- setClass("mcmcoutputpermfixpost",
                                   contains = c("mcmcpermfixpost", 
                                                "mcmcoutputfixpost"),
                                   validity = function(object) 
                                   {
                                       ## else: OK
                                       TRUE
                                   }
)

setMethod("initialize", "mcmcoutputpermfixpost",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), logperm = list(), 
                   postperm = list()) 
          {
              .Object@M         <- mcmcoutput@M
              .Object@burnin    <- mcmcout@burnin
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@post      <- mcmcoutput@post
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object@logperm   <- logperm
              .Object@postperm  <- postperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixpost",
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
              cat("     post        : List of",
                  length(object@post), "\n")
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     postperm    : List of",
                  length(object@postperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod("plot", signature(x = "mcmcoutputpermfixpost", 
                            y = "ANY"), 
          function(x, y = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permtraces.Poisson(x, y)
              }	
              ## log ##
              .permtraces.Log(x, y)
          }
)

setMethod("plotHist", signature(x = "mcmcoutputpermfixpost", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .permhist.Poisson(x, dev)
              }	
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputpermfixpost", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              if(x@model@dist == "poisson") {
                  .permdens.Poisson(x, dev)
              }	
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputpermfixpost",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              if (x@model@dist == "poisson") {
                  .permpointproc.Poisson(x, dev)
              }
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputpermfixpost",
                                   dev  = "ANY"),
          function(x, dev, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permsamprep.Poisson(x, dev)
              }
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputpermfixpost",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              if (x@model@dist == "poisson") {
                  .permpostdens.Poisson(x, dev)
              }
          }
)

