setClass("mcmcpermfix", 
         representation(
                        Mperm       = "integer",
                        parperm     = "list"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

## Define Classes in inheritance structure ##
setClass("mcmcoutputpermfix",
         contains = c("mcmcpermfix", "mcmcoutputfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermfix", 
          function(.Object, mcmcoutput, Mperm, parperm) {
              .Object@M         <- mcmcoutput@M
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfix",
          function(object) {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of ", 
                  length(object@par), "\n")
              cat("     log         : List of ", 
                  length(object@log), "\n")
              cat("     Mperm       : ", object@Mperm, "\n")
              cat("     parperm     : List of ", 
                  length(object@parperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

### --- Class 'mcmcoutputpermfixhier' --- ###
setClass("mcmcoutputpermfixhier",
         contains = c("mcmcpermfix", "mcmcoutputfixhier"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermfixhier",
          function(.Object, mcmcoutput) {
              .Object@M         <- mcmcoutput@M
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@hyper     <- mcmcoutput@hyper
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixhier",
          function(object) {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of ", 
                  length(object@par), "\n")
              cat("     log         : List of ", 
                  length(object@log), "\n")
              cat("     hyper       : List of ",
                  length(object@hyper), "\n")
              cat("     Mperm       : ", object@Mperm, "\n")
              cat("     parperm     : List of ", 
                  length(object@parperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

## Getters ##
setGeneric("getMPerm", function(object) standardGeneric("getMPerm"))
setMethod("getMPerm", "mcmcpermfix", 
          function(object) {
              w
              return(object@Mperm)
          }
)

setGeneric("getParperm", function(object) standardGeneric("getParperm"))
setMethod("getParperm", "mcmcpermfix", 
          function(object) {
              return(object@parperm)
          }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
