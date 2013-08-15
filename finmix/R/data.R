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
"data" <- function(y = matrix(), N, r, S = matrix(), 
                   bycolumn = TRUE, name = character(), 
			       type = "continuous", sim = FALSE, 
                   exp = matrix(), T = matrix()) {
		storage.mode(T) <- "integer"
		storage.mode(exp) <- "integer"
		
		has.data <- !all(is.na(y))
		if (missing(N) && has.data) {
			if(bycolumn) {
				N <- nrow(y)
			} else {
				N <- col(y)
			}
        } else if (missing(N) && !has.data) {
				N <- as.integer(1)
	    } else {
            N <- as.integer(N)
        }
		
		if (missing(r) && has.data) {
			if (bycolumn) {
                r <- ncol(y)
            } else {
                r <- nrow(y)
            } 
        } else if (missing(r) && !has.data){
            r <- as.integer(1)
        } else {
            r <- as.integer(r)
        }
		object <- new("data", y = y, N = N, r = r, S = S, 
                    bycolumn = bycolumn, name = name, type = type,
				    sim = sim, exp = exp, T = T)
		return(object)
}

setMethod("plot", "data", 
          function(x, y, ..., dev = TRUE){
              if (x@bycolumn) {
                  datam <- x@y
              } else {
                  datam <- t(x@y)
              }
              if (.check.grDevice() && dev) {
                  dev.new(title = "Histograms")
              }
              if (x@type == "discrete") {
                  ## univariate distributions ##
                  ## only univariates are implemented ##
                  barplot(table(datam), col = "gray65", 
                          border = "white", cex = 0.7,
                          cex.axis = 0.7, xlab = "", main = "",
                          cex.lab = 0.7)
                  if (!is.null(colnames(datam))) {
                      col.names <- colnames(datam)
                  } else {
                      col.names <- c("")
                  }
                  mtext(side = 1, col.names, cex = 0.7,
                        line = 3)     
              } else { ## continuous
                  if (x@r == 1) {
                      .symmetric.Hist(datam, colnames(datam))
                  } else if (x@r == 2) {        
                      .symmetric.Hist(datam, colnames(datam))
                      if (.check.grDevice() && dev) {
                          dev.new(title = "Contour plot")
                      }
                      par(mfrow = c(1, 2), mar = c(2, 2, 2, 3),
                          oma = c(4, 5, 1, 5))
                      plot(datam[, 1], datam[, 2], col = "gray47",
                           cex = 0.7, cex.axis = 0.7,
                           pch = 20, xlab = "", ylab = "", 
                           main = "")
                      mtext(side = 1, colnames(datam)[1],
                            cex = 0.7, line = 3)
                      mtext(side = 2, colnames(datam)[2],
                            cex = 0.7, line = 3)
                      d <- bkde2D(datam, 
                                  bandwidth = c(sd(datam[, 1]), 
                                                sd(datam[, 2])))
                      contour(d$x1, d$x2, d$fhat, col = "gray47",
                              cex = 0.7, cex.axis = 0.7, 
                              xlab = "", ylab = "")
                      mtext(side = 1, colnames(datam)[1],
                            cex = 0.7, line = 3)
                      mtext(side = 2, colnames(datam)[2],
                            cex = 0.7, line = 3)
                      if (.check.grDevice() && dev) {
                          dev.new("Perspective plot")
                      }
                      if (!is.null(colnames(datam))) {
                          col.names <- colnames(datam)
                      } else {
                          col.names <- c("", "")
                      }
                      persp(d$x1, d$x2, d$fhat, main = "",
                            xlab = col.names[1], ylab = col.names[2], 
                            zlab = "", col = "gray65", 
                            border = "gray47", theta = 55, phi = 30,
                            expand = 0.5, lphi = 190, ltheta = 90,
                            r = 40, d = 0.1, cex = 0.7, cex.axis = 0.7, 
                            cex.lab = 0.7, ticktype = "detailed")                   
                  } else { ## multivariate distribution
                      .symmetric.Hist(datam, colnames(datam))
                      if (.check.grDevice() && dev) {
                          dev.new(title = "Pairs")
                      }
                      pairs(datam, col = "gray47", pch = 20, 
                            cex = 0.7, cex.axis = 0.7, cex.labels = 1.3)                      
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
setMethod("getY", "data", function(object) {
				return(object@y)
			}
)
 
setMethod("getN", "data", function(object) {
				return(object@N)
			}
)

setMethod("getR", "data", function(object) {
				return(object@r)
			}
)

setMethod("getS", "data", function(object) {
				return(object@S)
			}
)

setMethod("getBycolumn", "data", function(object) {
					return(object@bycolumn)
				}
)

setMethod("getName", "data", function(object) {
					return(object@name)
				}
)

setMethod("getType", "data", function(object) {
					return(object@type)
				}
)

setMethod("getSim", "data", function(object) {
					return(object@sim)
				}
)

setMethod("getExp", "data", function(object) {
					return(object@exp)
				}
)

setMethod("getT", "data", function(object) {
					return(object@T)
				}
)

## Setters ##
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

setReplaceMethod("setN", "data", function(object, value) {
					object@N <- as.integer(value)
					validObject(object)
					return(object)
				}
)

setReplaceMethod("setR", "data", function(object, value) {
					object@r <- as.integer(value)
					validObject(object)
					return(object)
				}
)

setReplaceMethod("setS", "data", function(object, value) {
					object@S <- value
					validObject(object)
					return(object)
				}
)

setReplaceMethod("setBycolumn", signature(object = "data", value = "logical"), 
                 function(object, value) {						
						.Object@bycolumn <- value
						validObject(object)
						return(object)
					}
)  

setReplaceMethod("setName", "data", function(object, value) {
						object@name <- value
						validObject(object)
						return(object)
					}
)

setReplaceMethod("setType", "data", function(object, value) {
						object@type <- value
						validObject(object)
						return(object)
					}
)

setReplaceMethod("setSim", "data", function(object, value) {
						object@sim <- value
						validObject(object)
						return(object)
					}
)

setReplaceMethod("setExp", "data", function(object, value) {
						storage.mode(value) <- "integer"
						object@exp <- value
						validObject(object)
						return(object)
					}
)

setReplaceMethod("setT", "data", function(object, value) {
						storage.mode(value) <- "integer"
						object@T <- value
						validObject(object)
						return(object)
					}
)
