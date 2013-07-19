setClass("mcmcpermind", 
         representation(weightperm  = "array",
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
                   weightperm, STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
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
              cat("     par         : List of ", 
                  length(object@par), "\n")
              cat("     log         : List of ", 
                  length(object@log), "\n")
              cat("     ST          : ", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              cat("     S           : ", 
                  paste(dim(object@S), collapse = "x"), "\n")
              cat("     NK          : ",
                  paste(dim(object@NK), collapse = "x"), "\n")
              cat("     clust       : ",
                  paste(dim(object@clust), collapse = "x"), "\n")
              cat("     Mperm       : ", object@Mperm, "\n")
              cat("     parperm     : List of ", 
                  length(object@parperm), "\n")
              cat("     weightperm  : ",
                  paste(dim(object@weightperm), collapse = "x"), "\n")
              cat("     STperm      : ",
                  paste(dim(object@STperm), collapse = "x"), "\n")
              cat("     Sperm       : ",
                  paste(dim(object@Sperm), collapse = "x"), "\n")
              cat("     NKperm      : ", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)
setClass("mcmcoutputpermhier",
         contains = c("mcmcpermind", "mcmcoutputhier"),
         validity = function(object) {
             ## else: OK
             TRUE
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
setGeneric("getWeightPerm", function(object) standardGeneric("getWeightPerm"))
setMethod("getWeightPerm", "mcmcpermind", function(object) {
          return(object@weightperm)
})

## No setters as users are not intended to modify these ##
## obect.                                               ##
