.normalmodelmoments <- setClass("normalmodelmoments", 
                                representation(B   = "numeric",
                                               W   = "numeric",
                                               R   = "numeric"),
                                contains = c("cmodelmoments"),
                                validity = function(object) {
                                    ## else: OK
                                    TRUE
                                },
                                prototype(B = numeric(),
                                          W = numeric(),
                                          R = numeric()
                                          )
)

setMethod("initialize", "normalmodelmoments",
          function(.Object, ..., model) {
              .Object <- callNextMethod(.Object, ..., model = model)
              .Object <- generateMoments(.Object) 
              return(.Object)
          }
)

setMethod("generateMoments", "normalmodelmoments",
          function(object) 
          {
              .generateMomentsNormal(object)
          }
)

setMethod("show", "normalmodelmoments",
          function(object) {
              cat("Object 'modelmoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         :", 
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     higher      :",
                  paste(dim(object@higher), collapse = "x"), "\n")
              cat("     skewness    : Vector of", 
                  length(object@skewness), "\n")
              cat("     kurtosis    : Vector of", 
                  length(object@kurtosis), "\n")
              cat("     B           :", object@B, "\n")
              cat("     W           :", object@W, "\n")
              cat("     R           :", object@R, "\n")
              cat("     model       : Object of class",
                  class(object@model), "\n")
          }
)

## Getters ##
setMethod("getB", "normalmodelmoments", 
           function(object) 
           {
               return(object@B)
           }
)

setMethod("getW", "normalmodelmoments", 
           function(object) 
           {
               return(object@W)
           }
)

setMethod("getR", "normalmodelmoments", 
           function(object) 
           {
               return(object@R)
           }
)

## No setters as users are not intended to manipulate ##
## this object. ##

### Private functions
### This functions are not exported
".generateMomentsNormal" <- function(object) 
{
    mu              <- object@model@par$mu
    sigma           <- object@model@par$sigma
    weight          <- object@model@weight
    object@mean     <- sum(weight * mu)
    object@higher   <- .mixturemoments.normal(object@model, 
                                              4, object@mean)
    dimnames(object@higher) <- list(c("1st", "2nd", 
                                      "3rd", "4th"), "")
    object@var      <- array(object@higher[2], dim = c(1, 1))
    object@skewness <- object@higher[3]/object@higher[2]^1.5
    object@kurtosis <- object@higher[4]/object@higher[2]^2
    object@B        <- sum(weight * (mu - object@mean)^2)
    object@W        <- sum(weight * sigma)
    object@R        <- 1 - object@W/object@var[1]
    return(object)
}

