setClass("mcmcestfix",
         representation(dist        = "character",
                        K           = "integer",
                        indicmod    = "character",
                        map         = "list",
                        bml         = "list",
                        ieavg       = "list"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setMethod("show", "mcmcestfix", 
          function(object) {
              cat("Object 'mcmcest\n")
              cat("     K           :", object@K, "\n")
              cat("     indicmod    :", object@indicmod, 
                  "\n")
              cat("     map         : List of", 
                  length(object@map), "\n")
              cat("     bml         : List of",
                  length(object@bml), "\n")
              cat("     ieavg       : List of", 
                  length(object@ieavg), "\n")
          }
)

## Getters ##
## Generic set in 'model.R' ##
setMethod("getK", "mcmcestfix", 
          function(object) {
              return(object@K)
          }
)
## Generic set in 'model.R'
setMethod("getIndicmod", "mcmcestfix", 
          function(object) {
              return(object@indicmod)
          }
)

setGeneric("getMap", function(object) standardGeneric("getMap"))
setMethod("getMap", "mcmcestfix", 
           function(object) {
               return(object@map)
           }
)

setGeneric("getBml", function(object) standardGeneric("getBml"))
setMethod("getBml", "mcmcestfix",
          function(object) {
              return(object@bml)
          }
)

setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))
setMethod("getIeavg", "mcmcestfix", 
          function(object) {
              return(object@ieavg)
          }
)

## No setters as users are not intended to manipulate
## this object

