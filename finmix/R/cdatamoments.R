setClass("cdatamoments",
	representation(
		moments = "matrix",
		skewness = "matrix",
		kurtosis = "matrix",
		corr = "matrix",
		data = "data"),
	contains = "datamoments",
	prototype = list(
		moments = matrix(), 
		skewness = matrix(), 
		kurtosis = matrix(), 
		corr = matrix(), 
		data = data(), 
		mean = matrix(), 
		variance = matrix()),
	validity = function(.Object) {
		mom.moments <- .Object@moments
		mom.data <- .Object@data
		mom.r <- r(mom.data)
		if(mom.r != nrow(.Object@mean)) 
			return("Data dimension and dimension of the mean do not match.")
		if(mom.r != nrow(.Object@variance)) 
			return("Data dimension and dimension of the variance do not match.")
		if(mom.r != nrow(.Object@moments)) 
			return("Data dimension and dimension of the L=4 moments do not match.")
		if(mom.r != nrow(.Object@skewness)) 
			return("Data dimension and dimension of the skewness do not match.")
		if(mom.r != nrow(.Object@kurtosis)) 
			return("Data dimension and dimension of the kurtosis do not match.")
		## else: ok
		TRUE
	}
)


