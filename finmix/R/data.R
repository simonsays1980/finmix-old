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
                if (object@bycolumn) {
                    if (nrow(object@y) != object@N) {
                        warning("Slot 'N' does not correspond to 
                                the number of rows in slot 'y'.")
                    }
                    if (ncol(object@y) != object@r) {
                        warning("Slot 'r' does not correspond to 
                                the number of columns of slot 'y'.")
                    }
                } else {
                    if (ncol(object@y) != object@N) {
                        warning("Slot 'N' does not correspond to 
                                the number of columns of slot 'y'.")
                    }
                    if (nrow(object@y) != object@r) {
                        warning("Slot 'r' does not correspond to 
                                the number of rows of slot 'y'.")
                    }
                }
				if (has.T && any(object@T < 0)) 
					return("[Error] Repetitions 'T' smaller zero.")
				if (has.Exp && any(object@exp < 0)) 
					return("[Error] Exposures 'exp' smaller zero.")
				if (!(object@type %in% type.choices)) 
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
				object <- x
				lname <- length(object@name)
				if(object@type == "continuous") {
					if(object@bycolumn) {
						if(object@r > 1) {
							labels <- colnames(object@y)
							par(mfcol = c(ceiling(object@r/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:object@r) {
								xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
								hist(object@y[,i], xlab = xlab, main = "", ...) 
							}
							if( lname > 0) title(main = object@name, outer = TRUE, cex = 1.5)
							dev.new()
							if(r > 2) {
								par(mfcol=c(ceiling(choose(object@r,2)/2), 2), omi = c(0,0,0.3,0))
							}
							for(i in 1:(object@r - 1)) {
								for(j in (i + 1):object@r) {
									xlab <- ifelse(is.null(labels), paste("var ", i), 
											labels[i])
									ylab <- ifelse(is.null(labels), paste("var ", j), 
											labels[j]) 
									plot(object@y[,i], object@y[,j], xlab = xlab, 
										ylab = ylab, ...)
								}
							}
							if(r > 2) {
								if( lname > 0 ) title(main = object@name, outer = TRUE, cex = 1.5)
							}
							else{ ## r == 2
								title(main = object@name, outer = FALSE, cex = 1.5)
							}
							if(r > 2) {
								dev.new()
								main <- ifelse(lname == 0, "", object@name)
								pairs(object@y, main = main, ...)
							}
									 
						}
						else {
							main <- ifelse(lname == 0, "", object@name)
							xlab <- ifelse(is.null(colnames(object@y)), "",  
								colnames(object@y))
							hist(object@y, main = main , xlab = xlab, ...)
						}
					}
					else { ## by row 
						if(object@r > 1) {
                                                        labels <- rownames(object@y)
                                                        par(mfcol = c(ceiling(object@r/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:object@r) {
                                                                xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
                                                                hist(object@y[i,], xlab = xlab, main = "", ...)
                                                        }
                                                        if( lname > 0) title(main = object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        par(mfcol=c(ceiling(choose(object@r,2)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:(object@r - 1)) {
                                                                for(j in (i + 1):object@r) {
                                                                        xlab <- ifelse(is.null(labels), paste("var ", i),
                                                                                        labels[i])
                                                                        ylab <- ifelse(is.null(labels), paste("var ", j),
                                                                                        labels[j])
                                                                        plot(object@y[i,], object@y[j,], xlab = xlab,
                                                                                ylab = ylab, ...)
                                                                }
                                                        }
                                                        if( lname > 0 ) title(main = object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        main <- ifelse(lname == 0, "", object@name)
                                                        pairs(t(object)@y, main = main, ...)


                                                }
                                                else {
                                                        main <- ifelse(lname == 0, "", object@name)
                                                        xlab <- ifelse(is.null(rownames(object@y)), "",
                                                                colnames(object@y))
                                                        hist(t(object@y), main = main , xlab = xlab, ...)
						}
					}
							
					
				}
				else { # if type is 'discrete'
					if(object@bycolumn) {
						if(object@r > 1) {
							labels <- colnames(object@y)
							par(mfcol = c(ceiling(object@r/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:object@r) {
								xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
								barplot(table(object@y[,i]), xlab = xlab, main = "", ...) 
							}
							if( lname > 0) title(main = object@name, outer = TRUE, cex = 1.5)
							dev.new()
							par(mfcol=c(ceiling(choose(object@r,2)/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:(object@r - 1)) {
								for(j in (i + 1):ncol(object@y)) {
									xlab <- ifelse(is.null(labels), paste("var ", i), 
											labels[i])
									ylab <- ifelse(is.null(labels), paste("var ", j), 
											labels[j]) 
									plot(object@y[,i], object@y[,j], xlab = xlab, 
										ylab = ylab, ...)
								}
							}
							if( lname > 0 ) title(main = object@name, outer = TRUE, cex = 1.5)
							dev.new()
							main <- ifelse(lname == 0, "", object@name)
							pairs(object@y, main = main, ...)

									 
						}
						else { ## univariate data
							main <- ifelse(lname == 0, "", object@name)
							xlab <- ifelse(is.null(colnames(object@y)), "",  
								colnames(object@y))
							barplot(table(object@y), main = main , xlab = xlab, ...)
						}
					}
					else { ## ordered by row
						if(object@r> 1) {
                                                        labels <- rownames(object@y)
                                                        par(mfcol = c(ceiling(object@r/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:object@r) {
                                                                xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
                                                                barplot(table(t(object@y[i,])), xlab = xlab, main = "", ...)
                                                        }
                                                        if( lname > 0) title(main = object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        par(mfcol=c(ceiling(choose(object@r,2)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:(object@r - 1)) {
                                                                for(j in (i + 1):object@r) {
                                                                        xlab <- ifelse(is.null(labels), paste("var ", i),
                                                                                        labels[i])
                                                                        ylab <- ifelse(is.null(labels), paste("var ", j),
                                                                                        labels[j])
                                                                        plot(object@y[i,], object@y[j,], xlab = xlab,
                                                                                ylab = ylab, ...)
                                                                }
                                                        }
                                                        if( lname > 0 ) title(main = object@name, outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        main <- ifelse(lname == 0, "", object@r)
                                                        pairs(t(object@y), main = main, ...)


                                                }
                                                else { ## univariate discrete data ordered by row
                                                        main <- ifelse(lname == 0, "", object@name)
                                                        xlab <- ifelse(is.null(rownames(object@y)), "",
                                                                colnames(object@y))
                                                        barplot(table(t(object@y)), main = main , xlab = xlab, ...)
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
setGeneric("getY", function(object) standardGeneric("getY"))
setMethod("getY", "data", function(object) {
				return(object@y)
			}
)
 
setGeneric("getN", function(object) standardGeneric("getN"))
setMethod("getN", "data", function(object) {
				return(object@N)
			}
)

## Already set as generic in 'model.R' ##
setMethod("getR", "data", function(object) {
				return(object@r)
			}
)

setGeneric("getS", function(object) standardGeneric("getS"))
setMethod("getS", "data", function(object) {
				return(object@S)
			}
)

setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))
setMethod("getBycolumn", "data", function(object) {
					return(object@bycolumn)
				}
)

setGeneric("getName", function(object) standardGeneric("getName"))
setMethod("getName", "data", function(object) {
					return(object@name)
				}
)

setGeneric("getType", function(object) standardGeneric("getType"))
setMethod("getType", "data", function(object) {
					return(object@type)
				}
)

setGeneric("getSim", function(object) standardGeneric("getSim"))
setMethod("getSim", "data", function(object) {
					return(object@sim)
				}
)
setGeneric("getExp", function(object) standardGeneric("getExp"))
setMethod("getExp", "data", function(object) {
					return(object@exp)
				}
)
## Already set as Generic in 'model.R' ##
setMethod("getT", "data", function(object) {
					return(object@T)
				}
)
## Explicit usual R setter ##
setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))
setReplaceMethod("setY", "data", function(object, value) {
					if(object@bycolumn && NROW(value) == 1) {
						object@y <- t(value)
					}
					else {
						object@y <- value
					}
					if(object@bycolumn){
						object@N <- NROW(value)
						object@r <- NCOL(value)
					}
					else {
						object@N <- NCOL(value)
						object@r <- NROW(value)
					}
					validObject(object)
					return(object)
				}
)

setGeneric("setN<-", function(object, value) standardGeneric("setN<-"))
setReplaceMethod("setN", "data", function(object, value) {
					object@N <- as.integer(value)
					validObject(object)
					return(object)
				}
)

## Already set as generic in 'model.R' ##
setReplaceMethod("setR", "data", function(object, value) {
					object@r <- as.integer(value)
					validObject(object)
					return(object)
				}
)

setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))
setReplaceMethod("setS", "data", function(object, value) {
					object@S <- value
					validObject(object)
					return(object)
				}
)

setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))
setReplaceMethod("setBycolumn", signature(object = "data", value = "logical"), 
                 function(object, value) {						
						.Object@bycolumn <- value
						validObject(object)
						return(object)
					}
)  

setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))
setReplaceMethod("setName", "data", function(object, value) {
						object@name <- value
						validObject(object)
						return(object)
					}
)

setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))
setReplaceMethod("setType", "data", function(object, value) {
						object@type <- value
						validObject(object)
						return(object)
					}
)

setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))
setReplaceMethod("setSim", "data", function(object, value) {
						object@sim <- value
						validObject(object)
						return(object)
					}
)
setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))
setReplaceMethod("setExp", "data", function(object, value) {
						storage.mode(value) <- "integer"
						object@exp <- value
						validObject(object)
						return(object)
					}
)
## Already set as Generic in 'model.R' ##
setReplaceMethod("setT", "data", function(object, value) {
						storage.mode(value) <- "integer"
						object@T <- value
						validObject(object)
						return(object)
					}
)
