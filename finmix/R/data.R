## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

setClass("data", 
         representation (y           = "matrix",
		                 N           = "integer",
                   		 r           = "integer",
                		 S           = "matrix",
                		 bycolumn    = "logical",
                		 name        = "character",
		                 type        = "character",
		                 sim         = "logical",
		                 exp         = "matrix",
	    	             T           = "matrix"
	),	
	validity = function (object) {
				type.choices <- c("continuous", "discrete")
				has.T <- !all(is.na(object@T))
				has.Exp <- !all(is.na(object@exp))
				if(has.T && any(object@T < 0)) 
					return("[Error] Repetitions 'T' smaller zero.")
				if(has.Exp && any(object@exp < 0)) 
					return("[Error] Exposures 'exp' smaller zero.")
				if(!(object@type %in% type.choices)) 
					return("[Error] Unknown data 'type'. 'type' has to be 'continuous' or 'discrete'.")
				
				## else: ok
				TRUE
			}
)

## Constructor for the data class ##
"data" <- function(y. = matrix(), N., r., S. = matrix(), bycolumn. = TRUE, name. = character(), 
			type. = "continuous", sim. = FALSE, exp. = matrix(), T. = matrix()) {
		storage.mode(T.) <- "integer"
		storage.mode(exp.) <- "integer"
		
		has.data <- !all(is.na(y.))
		if (missing(N.) && has.data) {
			if(bycolumn.) {
				N. <- nrow(y.)
			} else {
				N. <- col(y.)
			}
        } else if (missing(N.) && !has.data) {
				N. <- as.integer(1)
	    } else {
            N. <- as.integer(N.)
        }
		
		if (missing(r.) && has.data) {
			if (bycolumn.) {
                r. <- ncol(y.)
            } else {
                r. <- nrow(y.)
            } 
        } else if (missing(r.) && !has.data){
            r. <- as.integer(1)
        } else {
            r. <- as.integer(r.)
        }
		data <- new("data", y = y., N = N., r = r., S = S., bycolumn = bycolumn., name = name., type = type.
				, sim = sim., exp = exp., T = T.)
		return(data)
}

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", "data", function(x, y, ..., deparse.level=1) {
				.Object <- x
				lname <- length(.Object@name)
				if(.Object@type == "continuous") {
					if(.Object@bycolumn) {
						if(.Object@r > 1) {
							labels <- colnames(.Object@y)
							par(mfcol = c(ceiling(.Object@r/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:.Object@r) {
								xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
								hist(.Object@y[,i], xlab = xlab, main = "", ...) 
							}
							if( lname > 0) title(main = .Object@name, outer = TRUE, cex = 1.5)
							dev.new()
							if(r > 2) {
								par(mfcol=c(ceiling(choose(.Object@r,2)/2), 2), omi = c(0,0,0.3,0))
							}
							for(i in 1:(.Object@r - 1)) {
								for(j in (i + 1):.Object@r) {
									xlab <- ifelse(is.null(labels), paste("var ", i), 
											labels[i])
									ylab <- ifelse(is.null(labels), paste("var ", j), 
											labels[j]) 
									plot(.Object@y[,i], .Object@y[,j], xlab = xlab, 
										ylab = ylab, ...)
								}
							}
							if(r > 2) {
								if( lname > 0 ) title(main = .Object@name, outer = TRUE, cex = 1.5)
							}
							else{ ## r == 2
								title(main = .Object@name, outer = FALSE, cex = 1.5)
							}
							if(r > 2) {
								dev.new()
								main <- ifelse(lname == 0, "", .Object@name)
								pairs(.Object@y, main = main, ...)
							}
									 
						}
						else {
							main <- ifelse(lname == 0, "", .Object@name)
							xlab <- ifelse(is.null(colnames(.Object@y)), "",  
								colnames(.Object@y))
							hist(.Object@y, main = main , xlab = xlab, ...)
						}
					}
					else { ## by row 
						if(.Object@r > 1) {
                                                        labels <- rownames(.Object@y)
                                                        par(mfcol = c(ceiling(.Object@r/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:.Object@r) {
                                                                xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
                                                                hist(.Object@y[i,], xlab = xlab, main = "", ...)
                                                        }
                                                        if( lname > 0) title(main = .Object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        par(mfcol=c(ceiling(choose(.Object@r,2)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:(.Object@r - 1)) {
                                                                for(j in (i + 1):.Object@r) {
                                                                        xlab <- ifelse(is.null(labels), paste("var ", i),
                                                                                        labels[i])
                                                                        ylab <- ifelse(is.null(labels), paste("var ", j),
                                                                                        labels[j])
                                                                        plot(.Object@y[i,], .Object@y[j,], xlab = xlab,
                                                                                ylab = ylab, ...)
                                                                }
                                                        }
                                                        if( lname > 0 ) title(main = .Object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        main <- ifelse(lname == 0, "", .Object@name)
                                                        pairs(t(.Object)@y, main = main, ...)


                                                }
                                                else {
                                                        main <- ifelse(lname == 0, "", .Object@name)
                                                        xlab <- ifelse(is.null(rownames(.Object@y)), "",
                                                                colnames(.Object@y))
                                                        hist(t(.Object@y), main = main , xlab = xlab, ...)
						}
					}
							
					
				}
				else { # if type is 'discrete'
					if(.Object@bycolumn) {
						if(.Object@r > 1) {
							labels <- colnames(.Object@y)
							par(mfcol = c(ceiling(.Object@r/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:.Object@r) {
								xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
								barplot(table(.Object@y[,i]), xlab = xlab, main = "", ...) 
							}
							if( lname > 0) title(main = .Object@name, outer = TRUE, cex = 1.5)
							dev.new()
							par(mfcol=c(ceiling(choose(.Object@r,2)/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:(.Object@r - 1)) {
								for(j in (i + 1):ncol(.Object@y)) {
									xlab <- ifelse(is.null(labels), paste("var ", i), 
											labels[i])
									ylab <- ifelse(is.null(labels), paste("var ", j), 
											labels[j]) 
									plot(.Object@y[,i], .Object@y[,j], xlab = xlab, 
										ylab = ylab, ...)
								}
							}
							if( lname > 0 ) title(main = .Object@name, outer = TRUE, cex = 1.5)
							dev.new()
							main <- ifelse(lname == 0, "", .Object@name)
							pairs(.Object@y, main = main, ...)

									 
						}
						else { ## univariate data
							main <- ifelse(lname == 0, "", .Object@name)
							xlab <- ifelse(is.null(colnames(.Object@y)), "",  
								colnames(.Object@y))
							barplot(table(.Object@y), main = main , xlab = xlab, ...)
						}
					}
					else { ## ordered by row
						if(.Object@r> 1) {
                                                        labels <- rownames(.Object@y)
                                                        par(mfcol = c(ceiling(.Object@r/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:.Object@r) {
                                                                xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
                                                                barplot(table(t(.Object@y[i,])), xlab = xlab, main = "", ...)
                                                        }
                                                        if( lname > 0) title(main = .Object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        par(mfcol=c(ceiling(choose(.Object@r,2)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:(.Object@r - 1)) {
                                                                for(j in (i + 1):.Object@r) {
                                                                        xlab <- ifelse(is.null(labels), paste("var ", i),
                                                                                        labels[i])
                                                                        ylab <- ifelse(is.null(labels), paste("var ", j),
                                                                                        labels[j])
                                                                        plot(.Object@y[i,], .Object@y[j,], xlab = xlab,
                                                                                ylab = ylab, ...)
                                                                }
                                                        }
                                                        if( lname > 0 ) title(main = .Object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        main <- ifelse(lname == 0, "", .Object@r)
                                                        pairs(t(.Object@y), main = main, ...)


                                                }
                                                else { ## univariate discrete data ordered by row
                                                        main <- ifelse(lname == 0, "", .Object@name)
                                                        xlab <- ifelse(is.null(rownames(.Object@y)), "",
                                                                colnames(.Object@y))
                                                        barplot(table(t(.Object@y)), main = main , xlab = xlab, ...)
						}
					}
				}
			}
)

setMethod("show", "data", 
          function(object) {

              has.S <- !all(is.na(object@S))
			  has.exp <- !all(is.na(object@exp))
			  has.T <- !all(is.na(object@T))
			  name <- ifelse(length(object@name) == 0, "data", object@name)
              cat("Object '", name, "'\n", sep = "")
              cat("     class       :", class(object), "\n")
              cat("     y           :", 
                  paste(dim(object@y), collapse = "x"), "\n")
              cat("     bycolumn    :", object@bycolumn, "\n")
              cat("     N           :", object@N, "\n")
              cat("     r           :", object@r, "\n")
              if (has.S) {
                  cat("     S           :", 
                      paste(dim(object@S), collapse = "x"), "\n")
              }
              cat("     type        :", object@type, "\n")
              cat("     sim         :", object@sim, "\n")
              if (has.exp) {
                  cat("     exp         :", 
                      paste(dim(object@exp), collapse = "x"), "\n")
              }
              if (has.T) {
                  cat("     T           :", 
                      paste(dim(object@T), collapse = "x"), "\n")
              }
          }
)
## Setters and Getters as a user interface to manipulate the slots ## 
## Combined Getter and Setter ##
setGeneric("getY", function(.Object) standardGeneric("getY"))
setMethod("getY", "data", function(.Object) {
				return(.Object@y)
			}
)
 
setGeneric("getN", function(.Object) standardGeneric("getN"))
setMethod("getN", "data", function(.Object) {
				return(.Object@N)
			}
)

## Already set as generic in 'model.R' ##
setMethod("getR", "data", function(object) {
				return(object@r)
			}
)

setGeneric("getS", function(.Object) standardGeneric("getS"))
setMethod("getS", "data", function(.Object) {
				return(.Object@S)
			}
)

setGeneric("getBycolumn", function(.Object) standardGeneric("getBycolumn"))
setMethod("getBycolumn", "data", function(.Object) {
					return(.Object@bycolumn)
				}
)

setGeneric("getName", function(.Object) standardGeneric("getName"))
setMethod("getName", "data", function(.Object) {
					return(.Object@name)
				}
)

setGeneric("getType", function(.Object) standardGeneric("getType"))
setMethod("getType", "data", function(.Object) {
					return(.Object@type)
				}
)

setGeneric("getSim", function(.Object) standardGeneric("getSim"))
setMethod("getSim", "data", function(.Object) {
					return(.Object@sim)
				}
)
setGeneric("getExp", function(.Object) standardGeneric("getExp"))
setMethod("getExp", "data", function(.Object) {
					return(.Object@exp)
				}
)
## Already set as Generic in 'model.R' ##
setMethod("getT", "data", function(.Object) {
					return(.Object@T)
				}
)
## Explicit usual R setter ##
setGeneric("setY<-", function(.Object, value) standardGeneric("setY<-"))
setReplaceMethod("setY", "data", function(.Object, value) {
					if(.Object@bycolumn && NROW(value) == 1) {
						.Object@y <- t(value)
					}
					else {
						.Object@y <- value
					}
					if(.Object@bycolumn){
						.Object@N <- NROW(value)
						.Object@r <- NCOL(value)
					}
					else {
						.Object@N <- NCOL(value)
						.Object@r <- NROW(value)
					}
					validObject(.Object)
					return(.Object)
				}
)

setGeneric("setN<-", function(.Object, value) standardGeneric("setN<-"))
setReplaceMethod("setN", "data", function(.Object, value) {
					.Object@N <- as.integer(value)
					validObject(.Object)
					return(.Object)
				}
)

## Already set as generic in 'model.R' ##
setReplaceMethod("setR", "data", function(.Object, value) {
					.Object@r <- as.integer(value)
					validObject(.Object)
					return(.Object)
				}
)

setGeneric("setS<-", function(.Object, value) standardGeneric("setS<-"))
setReplaceMethod("setS", "data", function(.Object, value) {
					.Object@S <- value
					validObject(.Object)
					return(.Object)
				}
)

setGeneric("setBycolumn<-", function(.Object, value) standardGeneric("setBycolumn<-"))
setReplaceMethod("setBycolumn", signature(.Object = "data", value = "logical"), 
                 function(.Object, value) {						
						.Object@bycolumn <- value
						validObject(.Object)
						return(.Object)
					}
)  

setGeneric("setName<-", function(.Object, value) standardGeneric("setName<-"))
setReplaceMethod("setName", "data", function(.Object, value) {
						.Object@name <- value
						validObject(.Object)
						return(.Object)
					}
)

setGeneric("setType<-", function(.Object, value) standardGeneric("setType<-"))
setReplaceMethod("setType", "data", function(.Object, value) {
						.Object@type <- value
						validObject(.Object)
						return(.Object)
					}
)

setGeneric("setSim<-", function(.Object, value) standardGeneric("setSim<-"))
setReplaceMethod("setSim", "data", function(.Object, value) {
						.Object@sim <- value
						validObject(.Object)
						return(.Object)
					}
)
setGeneric("setExp<-", function(.Object, value) standardGeneric("setExp<-"))
setReplaceMethod("setExp", "data", function(.Object, value) {
						storage.mode(value) <- "integer"
						.Object@exp <- value
						validObject(.Object)
						return(.Object)
					}
)
## Already set as Generic in 'model.R' ##
setReplaceMethod("setT", "data", function(.Object, value) {
						storage.mode(value) <- "integer"
						.Object@T <- value
						validObject(.Object)
						return(.Object)
					}
)
