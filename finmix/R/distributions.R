"dstud" <- function(x, mu, sigma, df) {
	
	fun <- gamma((df + 1)/2)/(gamma(df/2)*sqrt(df * pi * sigma)) * (1 + (x - mu)^2/(df * sigma))^(-(df + 1)/2)

	return(fun)
}
