.normultmodelmoments <- setClass("normultmodelmoments",
                                 representation(B   = "array",
                                                W   = "array",
                                                Rdet= "numeric",
                                                Rtr = "numeric",
                                                corr= "array"),
                                 contains = c("cmodelmoments"),
                                 validity = function(object) {
                                     ## else: OK
                                     TRUE
                                 },
                                 prototype(
                                           B    = array(),
                                           W    = array(),
                                           Rdet = numeric(),
                                           Rtr  = numeric(),
                                           corr = array()
                                           )
)

setMethod("initialize", "normultmodelmoments",
          function(.Object, ..., model) {
              .Object <- callNextMethod(.Object,..., model = model)
              generateMoments(.Object)
          }
)

setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))
setMethod("generateMoments", "normultmodelmoments",
          function(object) {
              .generateMomentsNormult(object)
          }
)

setMethod("show", "normultmodelmoments", 
          function(object) {
              cat("Object 'modelmoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         :",
                  paste(dim(object@var), collapse = "x"),
                  "\n")
              cat("     higher      :", 
                  paste(dim(object@higher), collapse = "x"),
                  "\n")
              cat("     skewness    : Vector of", 
                  length(object@skewness), "\n")
              cat("     kurtosis    : Vector of",
                  length(object@kurtosis), "\n")
              cat("     B           :",
                  paste(dim(object@B), collapse = "x"), "\n")
              cat("     W           :",
                  paste(dim(object@W), collapse = "x"), "\n")
              cat("     Rdet        :", object@Rdet, "\n")
              cat("     Rtr         :", object@Rtr, "\n")
              cat("     corr        :", 
                  paste(dim(object@corr), collapse = "x"),
                  "\n")
              cat("     model       : Object of class",
                  class(object@model), "\n")
          }
)

## Getters ##
setGeneric("getB", function(object) standardGeneric("getB"))
setMethod("getB", "normultmodelmoments", 
          function(object) {
              return(object@B)
          }
)

setGeneric("getW", function(object) standardGeneric("getW"))
setMethod("getW", "normultmodelmoments", 
          function(object) {
              return(object@W)
          }
)

setGeneric("getRdet", function(object) standardGeneric("getRdet"))
setMethod("getRdet", "normultmodelmoments", 
          function(object) {
              return(object@B)
          }
)

setGeneric("getRtr", function(object) standardGeneric("getRtr"))
setMethod("getRtr", "normultmodelmoments", 
          function(object) {
              return(object@B)
          }
)

setGeneric("getCorr", function(object) standardGeneric("getCorr"))
setMethod("getCorr", "normultmodelmoments",
          function(object) {
              return(object@corr)
          }
)

## No setters as users are not intended to manipulate ##
## this object ##

### private functions 
### these function are not exported 
".generateMomentsNormult" <- function(object) {
    mu          <- object@model@par$mu
    sigma       <- object@model@par$sigma
    weight      <- object@model@weight
    names       <- rep("", object@model@r)
    for (i in seq(1, object@model@r)) {
        names[i] <- paste("r=", i, sep = "")
    }
    object@mean <- apply(apply(mu, 1, '*', weight), 
                         2, sum, na.rm = TRUE)
    object@W    <- apply(sweep(sigma, MARGIN = 1, weight, '*'),
                         c(1,2), sum, na.rm = TRUE)
    object@var  <- object@W + apply(apply(mu, 2, tcrossprod, mu)
                                    , 1, '*', weight)
    object@var  <- object@var - object@mean %*% t(object@mean)
    diffm       <- mu - object@mean
    object@B    <- apply(apply(diffm, 1, tcrossprod, diffm),
                         1, '*', weight)
    cd          <- diag(1/diag(object@var)^.5)
    object@corr <- cd %*% object@var %*% cd
    object@Rtr  <- 1 - sum(diag(object@W))/sum(diag(object@var))
    object@Rdet <- 1 - det(object@W)/det(object@var)
    highm       <- array(0, dim = c(4, object@model@r))
    for(i in seq(1, object@model@r)) {
        marmodel    <- mixturemar(object@model, i)
        highm[, i]   <- t(.mixturemoments.normal(marmodel, 
                                                  4,
                                                  object@mean[i]))
    }
    names(object@mean)      <- names
    colnames(object@var)    <- names
    rownames(object@var)    <- names
    colnames(object@B)      <- names
    rownames(object@B)      <- names
    colnames(object@W)      <- names
    rownames(object@W)      <- names
    colnames(object@corr)   <- names
    rownames(object@corr)   <- names
    object@higher <- highm
    dimnames(object@higher) <- list(c("1st", "2nd", "3rd", "4th"), names)
    object@skewness <- object@higher[3, ]/object@higher[2,]^1.5
    object@kurtosis <- object@higher[4, ]/object@higher[2,]^2
    return(object)
}

