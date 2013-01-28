setClass("data", 
	representation (
		y = "matrix",
		r = "numeric",
		S = "matrix",
		bycolumn = "logical",
		name = "character",
		type = "character",
		sim = "logical",
		exp = "matrix",
		T = "matrix"
	),	
	validity = function (object) {
				type.choices <- c("continuous", "discrete")
				has.T <- !all(is.na(object@T))
				has.Exp <- !all(is.na(object@exp))
				if(has.T && any(object@T) < 0) 
					return("[Error] Repetitions 'T' smaller zero.")
				if(has.Exp && any(object@exp) < 0) 
					return("[Error] Exposures 'exp' smaller zero.")
				if(!(object@type %in% type.choices)) 
					return("[Error] Unknown data 'type'. 'type' has to be 'continuous' or 'discrete'.")
				
				## else: ok
				TRUE
			}
)

## Constructor for the data class ##
"data" <- function(y. = matrix(), r. = 1, S. = matrix(), bycolumn. = TRUE, name. = character(), 
			type. = "continuous", sim. = FALSE, exp. = matrix(), T. = matrix()) {
		storage.mode(T.) <- "integer"
		storage.mode(exp.) <- "integer"

		data <- new("data", y = y., r = r., S = S., bycolumn = bycolumn., name = name., type = type.
				, sim = sim., exp = exp., T = T.)
		
		return(data)
}

setMethod("plot", "data", function(x, ..., deparse.level=1) {
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

 setMethod("show", "data", function(object) {
					has.S <- !all(is.na(object@S))
					has.exp <- !all(is.na(object@exp))
					has.T <- !all(is.na(object@T))
					name <- ifelse(length(object@name) == 0, "", object@name)
					
					cat("Data object '", name, "'\n")
					cat("	Type		:", class(object), "\n")
					cat("	Data		:", paste(dim(object@y), collapse="x"), "\n")
					if(has.S) {
						classification <- paste(dim(object@S), collapse="x")
						cat("	Classifications	:", classification, "\n")
					}
					order <- ifelse(object@bycolumn, "by column", "by row")
					cat("	Order		:", order, "\n")
					cat("	Datatype	:", object@type, "\n")
					cat("	Simulated	:", object@sim, "\n")
					cat("	R (features)	:", object@r, "\n")
					if(has.exp) {
						exposures <- paste(dim(object@exp), collapse="x")
						cat("	Exposures	:", exposures, "\n")
					}
					if(has.T) {
						repetitions <- paste(dim(object@T), collapse="x")
						cat("	Repetitions	:", repetitions, "\n")
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
## Already set as generic in 'model.R' ##
setMethod("getR", "data", function(.Object) {
				return(.Object@r)
			}
)

setGeneric("getS", function(.Object) standardGeneric("getS"))
setMethod("getS", "data", function(.Object) {
				return(.Object@S)
			}
)

setGeneric("getByColumn", function(.Object) standardGeneric("getByColumn"))
setMethod("getByColumn", "data", function(.Object) {
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
					if(.Object@bycolumn && nrow(value) == 1) {
						.Object@y <- t(value)
					}
					else {
						.Object@y <- value
					}
					if(.Object@bycolumn){
						.Object@r <- ncol(value)
					}
					else {
						.Object@r <- nrow(value)
					}
					validObject(.Object)
					return(.Object)
				}
)
## Already set as generic in 'model.R' ##
setReplaceMethod("setR", "data", function(.Object, value) {
					.Object@r <- value
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

setGeneric("setByColumn<-", function(.Object, value) standardGeneric("setByColumn<-"))
setReplaceMethod("setByColumn", "data", function(.Object, value) {						
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
