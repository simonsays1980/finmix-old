#' @include data.R
#' @include groupmoments.R
setClass("sdatamoments",
	representation(
		group.moments = "groupmoments",
		WK = "array",
		var = "matrix",
		B = "matrix",
		W = "matrix",
		T = "matrix",
		Rtr = "numeric",
		Rdet = "numeric",
		s.data = "data"
		),
	validity = function(object) {
		mom.moments <- object@higher.moments
		mom.s.data <- object@s.data
		mom.r <- getR(mom.s.data)
		if(mom.r != ncol(object@mean)) 
			return("Data dimension and dimension of the mean do not match.")
		if(mom.r != ncol(object@variance)) 
			return("Data dimension and dimension of the variance do not match.")
		if(mom.r != ncol(mom.moments) || nrow(mom.moments) != 4) 
			return("Data dimension and dimension of the L=4 moments do not match.")
		if(mom.r != ncol(object@skewness)) 
			return("Data dimension and dimension of the skewness do not match.")
		if(mom.r != ncol(object@kurtosis)) 
			return("Data dimension and dimension of the kurtosis do not match.")
		if(!is.na(object@corr) && mom.r == 1)
			return("Data has dimension 1. No correlation matrix available.")
		if(!is.na(object@corr) && (mom.r != ncol(object@corr) || mom.r != nrow(object@corr)))
			return("Data dimension and dimension of correlation matrix do not match.")
		if(!is.na(object@kurtosis) && any(object@kurtosis < 0))
			return("Kurtosis is negative.")
		if(!is.na(object@variance) && any(object@variance < 0))
			return("Variance is negative")
		if(!is.na(mom.moments) && any(mom.moments[2,] < 0))
			return("Second higher moment is negative.")
		if(!is.na(mom.moments) && any(mom.moments[4,] < 0))
			return("Fourth higher moment is negative.")
		## else: ok
		TRUE
	}
)


setMethod("initialize", "sdatamoments", function(.Object, ..., data) {
						.Object@s.data <- data
						.Object@group.moments <- new("groupmoments", ..., data = data)
						smoments(.Object) <- data
						callNextMethod(.Object, ...)
						return(.Object)
					}
)

"smoments<-" <-  function(.Object, value) {
			              
				## Calculate moments ##
				## work only with data ordered by column ##
				if(!getByColumn(value)) {
					datam <- t(getY(value))
					classm <- t(getS(value))
				}
				else {
					datam <- getY(value)
					classm <- getS(value)
				}
				## Calculate within-group variability ##
				## 'WK' is an 1 x K array for r == 1 ##
				## 'WK' is an r x r x K array for r > 1 ##
				## Calculate the variance ##
				## 'var' is an 1 x K array for r == 1 ##
				## 'var' is an r x r x K array for r > 1 ##
				level.set <- as.numeric(levels(factor(classm)))
				K <- length(level.set)
				kmeans <- getGroupMean(getGroupMoments(.Object))
				r <- getR(getS.Data(.Object))
				nkm <- getNK(getGroupMoments(.Object))
				if(r == 1) {
					wkm <- array(0, dim = c(1, K))
					varm <- array(0, dim = c(1, K))
					for(i in 1:K) {
						group <- matrix(datam[which(classm == level.set[i]),])
						wkm[i] <- sum((group - kmeans[i])^2, na.rm = TRUE)
						varm[i] <- wkm[i]/kmeans[i]
					}
					.Object@WK <- wkm
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
					.Object@WK <- wkm
				}
				## Calculate the between-group variance ##
				## 'B' is an 1 x 1 matrix for r == 1 ##
				## 'B' is an r x r matrix for r > 1 ##
				bm <- matrix(0, ncol = r, nrow = 1)
				means <- as.matrix(apply(datam, 2, mean, na.rm = TRUE))
				out.prod <- matrix(0, ncol = r, nrow = r)
				for(i in 1:K) {
					out.prod <- out.prod + nkm[i] * (kmeans[,i] - means) %*% t(kmeans[,i] - means)  
				}	
				.Object@B <- as.matrix(out.prod)
				
				## Calculate the within-group heterogeneity ##
				## 'W' is an 1 x 1 matrix for r == 1 ##
				## 'W' is an r x r matrix for r > 1 ##
				if(r == 1){
					wm <- matrix(0, ncol = 1, nrow = 1)
					wm <- apply(getWK(.Object), 1, sum, na.rm = TRUE)
					.Object@W <- as.matrix(wm)
				}		
				else {
					wm <- matrix(0, ncol = r, nrow = r)
					for(i in 1:K) {
						wm <- wm + getWK(.Object)[,,i]
					}
					.Object@W <- wm
				}
				
				## Calculate total variance ##
				## 'T' is an 1 x 1 matrix for r == 1 ##
				## 'T' is an r x r matrix for r > 1 ##
				.Object@T <- as.matrix(getB(.Object) + getW(.Object))

				## Calculate coefficient of determination ##
				## 'Rtr' is an 1 x 1 numeric ##
				## 'Rdet' is an 1 x 1 numeric ## 
				w <- getW(.Object)
				t <- getT(.Object)
				.Object@Rtr <- 1 - (sum(diag(w), na.rm = TRUE) / sum(diag(t), na.rm = TRUE)) 
				.Object@Rdet <- 1 - (det(w)/det(t))
				
				return(.Object)
}

setMethod("show", "sdatamoments", function(object) {
						name.data <- getName(getS.Data(object))
						oname <- ifelse(length(name.data) == 0, "", name.data)
						cat("Moments object '", oname, "'\n")
						colnames <- colnames(getY(getS.Data(object)))
						gmoments <- getGroupMoments(object)
						cat("	group.moments:	  NK		: [", format(getNK(gmoments), trim = TRUE), "]\n")
						cat("			  group.mean	: [", format(getGroupMean(gmoments)[1,], trim = TRUE), "]\n")
						r <- getR(getS.Data(object))
						if(r > 1) {
							for(i in 2:r) {
								cat("			  	  [", format(getGroupMean(gmoments)[i,], trim = TRUE), "]\n")
							}
						}
						if(r == 1) {
							wk.array <- getWK(object)
							colnames(wk.array) <- c("K = 1", "K = 2")
							cat("	WK		: [", format(wk.array[1,], trim = TRUE), "]\n")
						}
						else {
							cat("	WK		: K = 1", "\n")
							for(i in 1:dim(getWK(object))[3]) {
								if(i > 1) {
									cat("		  	  K =", i,"\n")
								}
								mat <- getWK(object)[,,i]
								for(j in 1:nrow(mat)) {
									cat("		  	  	  [", format(mat[j,], trim = TRUE), "]\n")
								}
							}
						}
						cat("	B		: [", format(getB(object)[1,], trim = TRUE), "]\n")
						if(r > 1) {
							for(i in 2:r) {
								cat("		 	  [", format(getB(object)[i,], trim = TRUE), "]\n")
							}
						}
						cat("	W		: [", format(getW(object)[1,], trim = TRUE), "]\n")
						if(r > 1) {
							for(i in 2:r) {
								cat("		 	  [", format(getW(object)[i,], trim = TRUE), "]\n")
							}
						}
					        cat("	T		: [", format(getT(object)[1,], trim = TRUE), "]\n")
						if(r > 1) {
							for(i in 2:r) {
								cat("		 	  [", format(getT(object)[i,], trim = TRUE), "]\n")
							}
						}	
					 	cat("	Rtr		: [", format(getRtr(object), trim = TRUE), "]\n")
						cat("	Rdet		: [", format(getRdet(object), trim = TRUE), "]\n")	
					}
)

## Getters ##
setGeneric("getGroupMoments", function(.Object) standardGeneric("getGroupMoments"))
setMethod("getGroupMoments", "sdatamoments", function(.Object) {
						return(.Object@group.moments)
					}
)
setGeneric("getWK", function(.Object) standardGeneric("getWK"))
setMethod("getWK", "sdatamoments", function(.Object) {
						return(.Object@WK)
					}
)
setGeneric("getVar", function(.Object) standardGeneric("getVar"))
setMethod("getVar", "sdatamoments", function(.Object) {
				return(.Object@var)
			}
)
setGeneric("getB", function(.Object) standardGeneric("getB"))
setMethod("getB", "sdatamoments", function(.Object) {
						return(.Object@B)
					}
)
setGeneric("getW", function(.Object) standardGeneric("getW"))
setMethod("getW", "sdatamoments", function(.Object) {
						return(.Object@W)
					}
)
setGeneric("getT", function(.Object) standardGeneric("getT"))
setMethod("getT", "sdatamoments", function(.Object) {
						return(.Object@T)
					}
)
setGeneric("getRtr", function(.Object) standardGeneric("getRtr"))
setMethod("getRtr", "sdatamoments", function(.Object) {
						return(.Object@Rtr)
					}
)
setGeneric("getRdet", function(.Object) standardGeneric("getRdet"))
setMethod("getRdet", "sdatamoments", function(.Object) {
						return(.Object@Rdet)
					}
)
setGeneric("getS.Data", function(.Object) standardGeneric("getS.Data"))
setMethod("getS.Data", "sdatamoments", function(.Object) {
						return(.Object@s.data)				
					}
)
## Setters ##
setGeneric("setGroupMoments<-", function(.Object, value) standardGeneric("setGroupMoments<-"))
setReplaceMethod("setGroupMoments", "sdatamoments", function(.Object, value) {
							.Object@group.moments <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setWK<-", function(.Object, value) standardGeneric("setWK<-"))
setReplaceMethod("setWK", "sdatamoments", function(.Object, value) {
							.Object@WK <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setVar<-", function(.Object, value) standardGeneric("setVar<-"))
setReplaceMethod("setVar", "sdatamoments", function(.Object, value) {
							.Object@var <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setB<-", function(.Object, value) standardGeneric("setB<-"))			
setReplaceMethod("setB", "sdatamoments", function(.Object, value) {
							.Object@B <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setW<-", function(.Object, value) standardGeneric("setW<-"))			
setReplaceMethod("setW", "sdatamoments", function(.Object, value) {
							.Object@W <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setT<-", function(.Object, value) standardGeneric("setT<-"))			
setReplaceMethod("setT", "sdatamoments", function(.Object, value) {
							.Object@T <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setRtr<-", function(.Object, value) standardGeneric("setRtr<-"))			
setReplaceMethod("setRtr", "sdatamoments", function(.Object, value) {
							.Object@Rtr <- value
							validObject(.Object)
							return(.Object)
						}
)
setGeneric("setRdet<-", function(.Object, value) standardGeneric("setRdet<-"))
setReplaceMethod("setRdet", "sdatamoments", function(.Object, value) {
							.Object@Rdet <- value
							validObject(.Object)
							return(.Object)			
						}
)
setGeneric("setS.Data<-", function(.Object, value) standardGeneric("setS.Data<-"))
setReplaceMethod("setS.Data", "sdatamoments", function(.Object, value) {
							.Object@s.data <- value
							validObject(.Object)
							return(.Object)
						}
)
