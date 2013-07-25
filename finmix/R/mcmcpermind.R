setClass("mcmcpermind", 
         representation(weightperm  = "array",
                        entropyperm = "array",
                        STperm      = "array",
                        Sperm       = "array",
                        NKperm      = "array"
                        ),
         contains = c("mcmcpermfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

## Define Classes in inheritance structure ##
setClass("mcmcoutputpermbase",
         contains = c("mcmcpermind", "mcmcoutputbase"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermbase",
          function(.Object, mcmcoutput, Mperm, parperm,  
                   weightperm, logperm, entropyperm,
                   STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@weightperm    <- weightperm
              .Object@logperm       <- logperm
              .Object@entropyperm   <- entropyperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermbase", 
          function(object){
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              cat("     S           :", 
                  paste(dim(object@S), collapse = "x"), "\n")
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
              cat("     Sperm       :",
                  paste(dim(object@Sperm), collapse = "x"), "\n")
              cat("     NKperm      :", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

### --- Class 'mcmcoutputpermhier' --- ###
setClass("mcmcoutputpermhier",
         contains = c("mcmcpermind", "mcmcoutputhier"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermhier",
          function(.Object, mcmcoutput, Mperm, parperm,
                   weightperm, logperm, entropyperm, 
                   STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
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
          function(object){
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              cat("     S           :", 
                  paste(dim(object@S), collapse = "x"), "\n")
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
              cat("     Sperm       :",
                  paste(dim(object@Sperm), collapse = "x"), "\n")
              cat("     NKperm      :", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setClass("mcmcoutputpermpost",
         contains = c("mcmcpermind", "mcmcoutputpost"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setClass("mcmcoutputpermhierpost",
         contains = c("mcmcpermind", "mcmcoutputhier", "mcmcoutputpost"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setClassUnion("mcmcoutputperm",
              c("mcmcoutputpermfix",
                "mcmcoutputpermfixhier",
                "mcmcoutputpermfixpost",
                "mcmcoutputpermfixhierpost",
                "mcmcoutputpermbase",
                "mcmcoutputpermhier",
                "mcmcoutputpermpost",
                "mcmcoutputpermhierpost")
)

## Getters ##
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))
setMethod("getWeightperm", "mcmcpermind", function(object) {
          return(object@weightperm)
})

setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))
setMethod("getEntropyperm", "mcmcpermind", 
          function(object) {
              return(object@entropyperm)
          }
)

setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))
setMethod("getSTperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

setGeneric("getSperm", function(object) standardGeneric("getSperm"))
setMethod("getSperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))
setMethod("getNKperm", "mcmcpermind", 
          function(object) {
              return(object@STperm)
          }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
