setClass("mcmcpermfixpost",
         representation(postperm = "list"),
         contains = c("mcmcpermfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setClass("mcmcoutputpermfixpost",
         contains = c("mcmcpermfixpost", "mcmcoutputfixpost"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermfixpost",
          function(.Object, mcmcoutput) {
              .Object@M         <- mcmcoutput@M
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@post      <- mcmcoutput@post
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixpost",
          function(object) {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of ", 
                  length(object@par), "\n")
              cat("     log         : List of ", 
                  length(object@log), "\n")
              cat("     post        : List of ",
                  length(object@post), "\n")
              cat("     Mperm       : ", object@Mperm, "\n")
              cat("     parperm     : List of ", 
                  length(object@parperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

### --- Class 'mcmcoutputpermfixhierpost' --- ###
setClass("mcmcoutputpermfixhierpost",
         contains = c("mcmcpermfixpost", "mcmcoutputfixhierpost"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("initialize", "mcmcoutputpermfixhierpost",
          function(.Object, mcmcoutput) {
              .Object@M         <- mcmcoutput@M
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@hyper     <- mcmcoutput@hyper
              .Object@post      <- mcmcoutput@post
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixhierpost",
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
              cat("     post        : List of ",
                  length(object@post), "\n")
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
setGeneric("getPostperm", function(object) standardGeneric("getPostperm"))
setMethod("getPostperm", "mcmcpermfixpost", 
          function(object) {
              return(object@postperm)
          }
)

## No setters implemented as users are not intended to
## manipulate this object 
