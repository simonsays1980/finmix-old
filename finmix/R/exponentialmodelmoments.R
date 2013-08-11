.exponentialmodelmoments <- setClass("exponentialmodelmoments",
                                     representation(B   = "numeric",
                                                    W   = "numeric",
                                                    R   = "numeric"),
                                     contains = c("cmodelmoments"),
                                     validity = function(object) {
                                         ## else: OK
                                         TRUE
                                     },
                                     prototype(B    = numeric(),
                                               W    = numeric(),
                                               R    = numeric()
                                               )
)

setMethod("initialize", "exponentialmodelmoments",
          function(.Object, ..., model) 
          {
              .Object <- callNextMethod(.Object, ..., model = model)
              generateMoments(.Object)
          }
)

setMethod("generateMoments", "exponentialmodelmoments",
          function(object) 
          {
              .generateMomentsExponential(object)
          }
)

setMethod("show", "exponentialmodelmoments",
          function(object) {
              cat("Object 'modelmoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         :", 
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     higher      :", 
                  paste(dim(object@higher), collapse = "x"),
                  "\n")
              cat("     skewness    :", object@skewness, "\n")
              cat("     kurtosis    :", object@kurtosis, "\n")
              cat("     B           :", object@B, "\n")
              cat("     W           :", object@W, "\n")
              cat("     R           :", object@R, "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
          }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions 
### These functions are not exported 
".generateMomentsExponential" <- function(object) 
{
    lambda          <- object@model@par$lambda
    weight          <- object@model@weight
    object@mean     <- sum(weight * 1/lambda)
    highm <- .mixturemoments.exponential(object@model, 4, object@mean)
    dimnames(highm) <- list(c("1st", "2nd", "3rd", "4th"), "")
    object@higher   <- highm
    object@var      <- array(object@higher[2], dim = c(1, 1))
    object@skewness <- object@higher[3]/object@higher[2]^1.5
    object@kurtosis <- object@higher[4]/object@higher[2]^2
    object@W        <- sum(weight * 1/lambda^2)
    object@B        <- sum(weight * (1/lambda - object@mean)^2)
    object@R        <- 1 - object@W/object@var[1]
    return(object)
}

