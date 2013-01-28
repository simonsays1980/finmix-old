setClass("csdatamoments",
	representation(
		B = "matrix",
		W = "matrix",
		T = "matrix",
		Rtr = "numeric",
		Rdet = "numeric",
		data = "data"
		),
	contains = c("sdatamoments"),
	validity = function(object) {
		## else: ok
		TRUE
	}
)

".csdatamoments" <-  function(data) {
			              
				groupmoments <- groupmoments(data)
				## Calculate moments ##
				## work only with data ordered by column ##
				if(!data@bycolumn) {
					datam <- t(data@y)
					classm <- t(data@S)
				}
				else {
					datam <- data@y
					classm <- data@S
				}
				has.colnames <- !is.null(colnames(datam))
				## Calculate within-group variability ##
				## 'WK' is an 1 x K array for r == 1 ##
				## 'WK' is an r x r x K array for r > 1 ##
				## Calculate the variance ##
				## 'var' is an 1 x K array for r == 1 ##
				## 'var' is an r x r x K array for r > 1 ##
				level.set <- as.numeric(levels(factor(classm)))
				K <- length(level.set)
				kmeans <- groupmoments@group.mean
				r <- ncol(datam)
				nkm <- groupmoments@NK
				if(r == 1) {
					wkm <- array(0, dim = c(1, K))
					varm <- array(0, dim = c(1, K))
					for(i in 1:K) {
						group <- matrix(datam[which(classm == level.set[i]),])
						wkm[i] <- sum((group - kmeans[i])^2, na.rm = TRUE)
						varm[i] <- wkm[i]/kmeans[i]
					}
					dimnames(wkm) <- list(NULL, paste("k=", 1:K, sep=""))
					dimnames(varm) <- list(NULL, paste("k=", 1:K, sep=""))
					if(has.colnames) {
						rownames(wkm) <- colnames(datam)
						rownames(varm) <- colnames(datam)
					}
				}
				else { ## r > 1
					wkm <- array(0, c(r, r, K))
					varm <- array(0, c(r, r, K))
					for(i in 1:K) {
						kmeanm <- t(as.matrix(kmeans[,i])[, rep(1, nkm[i])])
						group <- as.matrix(datam[which(classm == level.set[i]),]) - kmeanm
						out.prod <- matrix(0, ncol  = r, nrow = r)
						for(j in 1:r) {
							out.prod <- out.prod + group[j,] %o% group[j,] 	
						}
						wkm[,,i] <- out.prod
						varm[,,i] <- out.prod/nkm[i] 
					}
					if(has.colnames) {
						dimnames(wkm) <- list(colnames(datam), colnames(datam), 
									paste("k=", 1:K, sep = ""))
						dimnames(varm) <- list(colnames(datam), colnames(datam), 
									paste("k=", 1:K, sep = ""))
					} 
					else {
						dimnames(wkm) <- list(NULL, NULL, paste("k=", 1:K, sep=""))
						dimnames(varm) <- list(NULL, NULL, paste("k=", 1:K, sep=""))
					}
				}
				## Calculate the between-group variance ##
				## 'B' is an 1 x 1 matrix for r == 1 ##
				## 'B' is an r x r matrix for r > 1 ##
				means <- as.matrix(apply(datam, 2, mean, na.rm = TRUE))
				out.prod <- matrix(0, ncol = r, nrow = r)
				for(i in 1:K) {
					out.prod <- out.prod + nkm[i] * (kmeans[,i] - means) %*% t(kmeans[,i] - means)  
				}	
				if(has.colnames) {
					colnames(out.prod) <- colnames(datam)
					rownames(out.prod) <- colnames(datam)
				}
				## Calculate the within-group heterogeneity ##
				## 'W' is an 1 x 1 matrix for r == 1 ##
				## 'W' is an r x r matrix for r > 1 ##
				if(r == 1){
					wm <- matrix(0, ncol = 1, nrow = 1)
					wm <- apply(wkm, 1, sum, na.rm = TRUE)
				}		
				else {
					wm <- matrix(0, ncol = r, nrow = r)
					for(i in 1:K) {
						wm <- wm + wkm[,,i]
					}
				}
				wm <- as.matrix(wm)
				if(has.colnames) {
                                        colnames(wm) <- colnames(datam)
                                        rownames(wm) <- colnames(datam)
                                }
				## Calculate total variance ##
				## 'T' is an 1 x 1 matrix for r == 1 ##
				## 'T' is an r x r matrix for r > 1 ##
				tm <- as.matrix(out.prod + wm)
				## Calculate coefficient of determination ##
				## 'Rtr' is an 1 x 1 numeric ##
				## 'Rdet' is an 1 x 1 numeric ##
				if(r > 1) { 
					rtr <- 1 - (sum(diag(wm), na.rm = TRUE) / sum(diag(tm), na.rm = TRUE)) 
					rdet <- 1 - (det(wm)/det(tm))
						
					csdatamoments <- new("csdatamoments", groupmoments = groupmoments, WK = wkm,
								var = varm, B = as.matrix(out.prod), W = as.matrix(wm),
								T = tm, Rdet = rdet, Rtr = rtr, data = data)
				}
					csdatamoments <- new("csdatamoments", groupmoments = groupmoments, WK = wkm,
                                                                var = varm, B = as.matrix(out.prod), W = as.matrix(wm),
                                                                T = tm, Rdet = numeric(), Rtr = numeric(), data = data)
				return(csdatamoments)
}

setMethod("show", "csdatamoments", function(object) {
						name.data <- getName(object@data)
						oname <- ifelse(length(name.data) == 0, "", name.data)
						cat("SDataMoments object '", oname, "'\n")
						gmoments <- object@groupmoments
						cat("	Type		:", class(object), "\n")
						cat("	Group Moments	:", class(object@groupmoments), "\n")
						cat("	WK		:", paste(dim(object@WK), collapse="x"), "\n")
						cat("	Var (within)	:", paste(dim(object@var), collapse="x"), "\n")
						cat("	B (var between)	:", paste(dim(object@B), collapse="x"), "\n")
						cat("	W (het within)	:", paste(dim(object@W), collapse="x"), "\n")
						cat("	T (var total)	:", paste(dim(object@T), collapse="x"), "\n")
						if(length(object@Rtr) > 0) {
						 	cat("	Rtr		: [", format(getRtr(object), trim = TRUE), "]\n")
						}
						if(length(object@Rdet) > 0) {
							cat("	Rdet		: [", format(getRdet(object), trim = TRUE), "]\n")	
						}
					}
)

## Getters ##
setMethod("getGroupMoments", "csdatamoments", function(.Object) {
						return(.Object@groupmoments)
					}
)
setMethod("getWK", "csdatamoments", function(.Object) {
						return(.Object@WK)
					}
)
setMethod("getVar", "csdatamoments", function(.Object) {
				return(.Object@var)
			}
)
setGeneric("getB", function(.Object) standardGeneric("getB"))
setMethod("getB", "csdatamoments", function(.Object) {
						return(.Object@B)
					}
)
setGeneric("getW", function(.Object) standardGeneric("getW"))
setMethod("getW", "csdatamoments", function(.Object) {
						return(.Object@W)
					}
)
## Already set as generic in 'model.R' ## 
setMethod("getT", "csdatamoments", function(.Object) {
						return(.Object@T)
					}
)
setGeneric("getRtr", function(.Object) standardGeneric("getRtr"))
setMethod("getRtr", "csdatamoments", function(.Object) {
						return(.Object@Rtr)
					}
)
setGeneric("getRdet", function(.Object) standardGeneric("getRdet"))
setMethod("getRdet", "csdatamoments", function(.Object) {
						return(.Object@Rdet)
					}
)
setMethod("getData", "csdatamoments", function(.Object) {
						return(.Object@data)				
					}
)

## Setters ##
## No setters, as it users are adviced not to manipulate moment objects ##
