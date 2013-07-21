setClass("mcmcestind",
         representation(eavg = "list"),
         contains = c("mcmcestfix"),
         validity = function(object) {
             ## else: OK
             TRUE
         }
)

setClassUnion("mcmcest", c("mcmcestfix",
                           "mcmcestind")
)

setMethod("show", "mcmcestind",
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
              cat("     eavg        : List of",
                  length(object@eavg), "\n")
          }
)

## Getters ##
setGeneric("getEavg", function(object) standardGeneric("getEavg"))
setMethod("getEavg", "mcmcestind", 
          function(object) {
              return(object@eavg)
          }
)

## No setters as users are not intended to manipulate 
## this object.
