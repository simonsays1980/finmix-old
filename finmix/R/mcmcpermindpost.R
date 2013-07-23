setClass("mcmcpermindpost",
         representation(postperm = "list"),
         contains = c("mcmcpermind"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)


## Class inheritance structure ##
setClass("mcmcoutputpermpost",
         contains = c("mcmcpermindpost", "mcmcoutputpost"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermpost",
          function(.Object, mcmcoutput, Mperm, parperm,
                   weightperm, postperm, STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@post          <- mcmcoutput@post
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@weightperm    <- weightperm
              .Object@postperm      <- postperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermpost", 
          function(object){
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     weight      :",
                  paste(dim(object@weight), collapse = "x"), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     post        : List of",
                  length(object@post), "\n")
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
              cat("     postperm    : List of",
                  length(object@postperm), "\n")
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

### --- Class 'mcmcoutputpermhierpost' --- ###
setClass("mcmcoutputpermhierpost",
         contains = c("mcmcpermindpost", "mcmcoutputhierpost"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermhierpost",
          function(.Object, mcmcoutput, Mperm, parperm,
                   weightperm, postperm, STperm, Sperm, NKperm) {
              .Object@M             <- mcmcoutput@M
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@hyper         <- mcmcoutput@hyper
              .Object@post          <- mcmcoutput@post
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@weightperm    <- weightperm
              .Object@postperm      <- postperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermhierpost", 
          function(object){
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     weight      :",
                  paste(dim(object@weight), collapse = "x"), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     post        : List of",
                  length(object@post), "\n")
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
              cat("     postperm    : List of",
                  length(object@postperm), "\n")
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
## Generic defined in 'mcmcpermfixpost.R' ##
setMethod("getPostperm", "mcmcpermindpost",
          function(object) {
              return(object@postperm)
          }
)

## No setters as users are not intended to manipulate 
## this objects
