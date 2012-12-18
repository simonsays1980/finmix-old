setClass("data", 
	representation (
		y = "matrix",
		r = "numeric",
		S = "matrix",
		bycolumn = "logical",
		name = "character",
		type = "character",
		sim = "logical",
		model = "model"),
	contains = "model",
	prototype = list(
			y = matrix(), 
			r = 1, 
			S = matrix(), 
			bycolumn = FALSE, 
			name = character(),
			type = "continuous",
			sim = FALSE,
			model = model()),
	validity = function (object) {
				data.y <- object@y
				data.r <- object@r
				if( !is.na(data.y) && ncol(data.y) == 1 && nrow(data.y) == 1 )
					return("Observations 'y' contain a single point. A finite mixture analysis cannot be conducted.")
				if( ncol(data.y) != data.r && nrow(data.y) != data.r )
					return("Dimensions of observations 'y' and dimension 'r' do not match.")
				data.S <- object@S
				if( !is.na(data.S) && (nrow(data.S) != nrow(data.y) || ncol(data.S) != ncol(data.y)) ) 
					return("Dimensions of classifications 'S' and observations 'y' do not match.")
				data.type <- object@type
				data.type.f <- strsplit(data.type, split = '')[[1]][1]
				if( data.type.f != "c" && data.type.f != "d" )
					return(paste("Data 'type' '", data.type, "' unknown. 'type' has to be either 'discrete' or 'continuous'."), sep="")
				## else: ok
				TRUE
			}
)

## Constructor for the data class ##
"data" <- function(y. = matrix(), r. = 1, S. = matrix(), bycolumn. = FALSE,
			name. = character(), type. = "continuous", sim. = FALSE, model. = model()) {
		data <- new("data", y = y., r = r., S = S., bycolumn = bycolumn., 
				name = name., type = type., sim = sim., model = model.)
		return(data)
}

## Initialization: Ensuring that type is either 'discrete' or 'continuous' ##
setMethod("initialize", "data", function(.Object, type., ...) {
					
					data.type.f <- strsplit(.Object@type, split = '')[[1]][1]	
					if( data.type.f == "d") .Object@type <- "discrete"
					else if( data.type.f == "c") .Object@type = "continuous"
					
					callNextMethod()			
				}
)

setMethod("plot", "data", function(x, ..., deparse.level=1) {
				.Object <- x
				lname <- length(name(.Object))
				if(type(.Object) == "continuous") {
					if(bycolumn(.Object)) {
						if(r(.Object) > 1) {
							labels <- colnames(y(.Object))
							par(mfcol = c(ceiling(r(.Object)/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:r(.Object)) {
								xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
								hist(y(.Object)[,i], xlab = xlab, main = "", ...) 
							}
							if( lname > 0) title(main = name(.Object), outer = TRUE, cex = 1.5)
							dev.new()
							par(mfcol=c(ceiling(choose(r(.Object),2)/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:(r(.Object) - 1)) {
								for(j in (i + 1):ncol(y(.Object))) {
									xlab <- ifelse(is.null(labels), paste("var ", i), 
											labels[i])
									ylab <- ifelse(is.null(labels), paste("var ", j), 
											labels[j]) 
									plot(y(.Object)[,i], y(.Object)[,j], xlab = xlab, 
										ylab = ylab, ...)
								}
							}
							if( lname > 0 ) title(main = name(.Object), outer = TRUE, cex = 1.5)
							dev.new()
							main <- ifelse(lname == 0, "", name(.Object))
							pairs(y(.Object), main = main, ...)

									 
						}
						else {
							main <- ifelse(lname == 0, "", name(.Object))
							xlab <- ifelse(is.null(colnames(y(.Object))), "",  
								colnames(y(.Object)))
							hist(y(.Object), main = main , xlab = xlab, ...)
						}
					}
					else {
						if(r(.Object) > 1) {
                                                        labels <- rownames(y(.Object))
                                                        par(mfcol = c(ceiling(r(.Object)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:r(.Object)) {
                                                                xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
                                                                hist(y(.Object)[i,], xlab = xlab, main = "", ...)
                                                        }
                                                        if( lname > 0) title(main = name(.Object), outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        par(mfcol=c(ceiling(choose(r(.Object),2)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:(r(.Object) - 1)) {
                                                                for(j in (i + 1):r(.Object)) {
                                                                        xlab <- ifelse(is.null(labels), paste("var ", i),
                                                                                        labels[i])
                                                                        ylab <- ifelse(is.null(labels), paste("var ", j),
                                                                                        labels[j])
                                                                        plot(y(.Object)[i,], y(.Object)[j,], xlab = xlab,
                                                                                ylab = ylab, ...)
                                                                }
                                                        }
                                                        if( lname > 0 ) title(main = name(.Object), outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        main <- ifelse(lname == 0, "", name(.Object))
                                                        pairs(t(y(.Object)), main = main, ...)


                                                }
                                                else {
                                                        main <- ifelse(lname == 0, "", name(.Object))
                                                        xlab <- ifelse(is.null(rownames(y(.Object))), "",
                                                                colnames(y(.Object)))
                                                        hist(t(y(.Object)), main = main , xlab = xlab, ...)
						}
					}
							
					
				}
				else { # if type is 'discrete'
					if(bycolumn(.Object)) {
						if(r(.Object) > 1) {
							labels <- colnames(y(.Object))
							par(mfcol = c(ceiling(r(.Object)/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:r(.Object)) {
								xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
								barplot(table(y(.Object)[,i]), xlab = xlab, main = "", ...) 
							}
							if( lname > 0) title(main = name(.Object), outer = TRUE, cex = 1.5)
							dev.new()
							par(mfcol=c(ceiling(choose(r(.Object),2)/2), 2), omi = c(0,0,0.3,0))
							for(i in 1:(r(.Object) - 1)) {
								for(j in (i + 1):ncol(y(.Object))) {
									xlab <- ifelse(is.null(labels), paste("var ", i), 
											labels[i])
									ylab <- ifelse(is.null(labels), paste("var ", j), 
											labels[j]) 
									plot(y(.Object)[,i], y(.Object)[,j], xlab = xlab, 
										ylab = ylab, ...)
								}
							}
							if( lname > 0 ) title(main = name(.Object), outer = TRUE, cex = 1.5)
							dev.new()
							main <- ifelse(lname == 0, "", name(.Object))
							pairs(y(.Object), main = main, ...)

									 
						}
						else {
							main <- ifelse(lname == 0, "", name(.Object))
							xlab <- ifelse(is.null(colnames(y(.Object))), "",  
								colnames(y(.Object)))
							barplot(table(y(.Object)), main = main , xlab = xlab, ...)
						}
					}
					else {
						if(r(.Object) > 1) {
                                                        labels <- rownames(y(.Object))
                                                        par(mfcol = c(ceiling(r(.Object)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:r(.Object)) {
                                                                xlab <- ifelse(is.null(labels), paste("var ", i), labels[i])
                                                                barplot(table(t(y(.Object)[i,])), xlab = xlab, main = "", ...)
                                                        }
                                                        if( lname > 0) title(main = name(.Object), outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        par(mfcol=c(ceiling(choose(r(.Object),2)/2), 2), omi = c(0,0,0.3,0))
                                                        for(i in 1:(r(.Object) - 1)) {
                                                                for(j in (i + 1):r(.Object)) {
                                                                        xlab <- ifelse(is.null(labels), paste("var ", i),
                                                                                        labels[i])
                                                                        ylab <- ifelse(is.null(labels), paste("var ", j),
                                                                                        labels[j])
                                                                        plot(y(.Object)[i,], y(.Object)[j,], xlab = xlab,
                                                                                ylab = ylab, ...)
                                                                }
                                                        }
                                                        if( lname > 0 ) title(main = name(.Object), outer = TRUE, cex = 1.5)
                                                        dev.new()
                                                        main <- ifelse(lname == 0, "", name(.Object))
                                                        pairs(t(y(.Object)), main = main, ...)


                                                }
                                                else {
                                                        main <- ifelse(lname == 0, "", name(.Object))
                                                        xlab <- ifelse(is.null(rownames(y(.Object))), "",
                                                                colnames(y(.Object)))
                                                        barplot(table(t(y(.Object))), main = main , xlab = xlab, ...)
						}
					}
				}
			}
)

setMethod("show", "data", function(object) {
					.Object <- object
					name <- ifelse(length(name(.Object)) == 0, "", name(.Object))
					cat("Data object '", name, "'\n")
					cat("	Type		:", class(.Object), "\n")
					cat("	Data		:", paste(dim(y(.Object)), collapse="x"), "\n")
					classification <- ifelse(is.na(S(.Object)), "", paste(dim(S(.Object)), collapse="x"))
					cat("	Classifications	:", classification, "\n")
					order <- ifelse(bycolumn(.Object), "by column", "by row")
					cat("	Order		:", order, "\n")
					cat("	Datatype	:", type(.Object), "\n")
					simulated <- ifelse(sim(.Object), "true", "false")
					cat("	Simulated	:", simulated, "\n")
					cat("	Model		:\n")
				}

)
## Setters and Getters as a user interface to manipulate the slots ## 
## Combined Getter and Setter ##
setGeneric("y", function(.Object) standardGeneric("y"))
setMethod("y", "data", function(.Object) {
				return(.Object@y)
			}
) 

setGeneric("r", function(.Object) standardGeneric("r"))
setMethod("r", "data", function(.Object) {
				return(.Object@r)
			}
)

setGeneric("S", function(.Object) standardGeneric("S"))
setMethod("S", "data", function(.Object) {
				return(.Object@S)
			}
)

setGeneric("bycolumn", function(.Object) standardGeneric("bycolumn"))
setMethod("bycolumn", "data", function(.Object) {
					return(.Object@bycolumn)
				}
)

setGeneric("name", function(.Object) standardGeneric("name"))
setMethod("name", "data", function(.Object) {
					return(.Object@name)
				}
)

setGeneric("type", function(.Object) standardGeneric("type"))
setMethod("type", "data", function(.Object) {
					return(.Object@type)
				}
)

setGeneric("sim", function(.Object) standardGeneric("sim"))
setMethod("sim", "data", function(.Object) {
					return(.Object@sim)
				}
)

setGeneric("data.model", function(.Object) standardGeneric("data.model"))
setMethod("data.model", "data", function(.Object) {
					return(.Object@model)
				}
)
## Explicit usual R setter ##
setGeneric("y<-", function(.Object, value) standardGeneric("y<-"))
setReplaceMethod("y", "data", function(.Object, value) {
					.Object@y <- value
					if( bycolumn(.Object) ) .Object@r <- ncol(value)
					else .Object@r <- nrow(value)
					validObject(.Object)
					return(.Object)
				}
)

setGeneric("r<-", function(.Object, value) standardGeneric("r<-"))
setReplaceMethod("r", "data", function(.Object, value) {
					.Object@r <- value
					validObject(.Object)
					return(.Object)
				}
)

setGeneric("S<-", function(.Object, value) standardGeneric("S<-"))
setReplaceMethod("S", "data", function(.Object, value) {
					.Object@S <- value
					validObject(.Object)
					return(.Object)
				}
)

setGeneric("bycolumn<-", function(.Object, value) standardGeneric("bycolumn<-"))
setReplaceMethod("bycolumn", "data", function(.Object, value) {
						bycolumn <- bycolumn(.Object)
						.Object@bycolumn <- value
						if( (value != bycolumn) && bycolumn) .Object@r <- nrow(y(.Object))
						else if((value != bycolumn) && !bycolumn) .Object@r <- ncol(y(.Object))
						validObject(.Object)
						return(.Object)
					}
)  

setGeneric("name<-", function(.Object, value) standardGeneric("name<-"))
setReplaceMethod("name", "data", function(.Object, value) {
						.Object@name <- value
						validObject(.Object)
						return(.Object)
					}
)

setGeneric("type<-", function(.Object, value) standardGeneric("type<-"))
setReplaceMethod("type", "data", function(.Object, value) {
						.Object@type <- value
						validObject(.Object)
						return(.Object)
					}
)

setGeneric("sim<-", function(.Object, value) standardGeneric("sim<-"))
setReplaceMethod("sim", "data", function(.Object, value) {
						.Object@sim <- value
						validObject(.Object)
						return(.Object)
					}
)

setGeneric("data.model<-", function(.Object, value) standardGeneric("data.model<-"))
setReplaceMethod("data.model", "data", function(.Object, value) {
						.Object@model <- value
						validObject(.Object)
						return(.Object)
					}
)

