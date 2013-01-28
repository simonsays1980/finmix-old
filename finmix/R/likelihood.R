"likelihood.normal" <- function(y, mu, sigma){
	N <- nrow(y)
	K <- ncol(mu)
	y <- matrix(y, nrow = N, ncol = K)
	
	err <- t(apply(y, 1, "-", mu))
	err <- t(apply(err^2, 1, "/", sigma))
	
	loglik <- -.5 * (log(2 * pi) + t(apply(err, 1, "+", log(sigma))))

	if(K == 1) {
		max.lik <- loglik
	}
	else { 
		max.lik <- apply(loglik, 1, max, na.rm = TRUE)
	}
	l.h <- exp(loglik -  max.lik)

	result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
	return(result)
}

"likelihood.student" <- function(y, mu, sigma, df) {
	N <- nrow(y)
	K <- ncol(mu)
	y <- matrix(y, nrow = N, ncol = K)

	err <- t(apply(y, 1, "-", mu))
	err <- t(apply(err^2, 1, "/", sigma))
	err <- t(apply(err, 1, "/", df))
	err <- t(apply(log(1 + err), 1, "*", (df + 1)/2))
	
	loglik <- lgamma((df + 1)/2) - lgamma(df/2) - .5 * (log(df * pi) + log(sigma))
	loglik <- t(apply(-err, 1, "+", loglik))

	if(K == 1) {
		max.lik <- loglik
	}
	else { 
		max.lik <- apply(loglik, 1, max, na.rm = TRUE)
	}
	l.h <- exp(loglik - max.lik)

	result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
	return(result)

}

"likelihood.exponential" <- function(y, lambda) {
	
	N <- nrow(y) 
	K <- ncol(lambda)
	y <- matrix(y, nrow = N, ncol = K)
	
	lambda <- matrix(lambda, nrow = N, ncol = K, byrow = TRUE)

	lambda <- apply(lambda, c(1,2), max, 0.0001)
	loglik <- log(lambda) - y * lambda
	max.lik <- apply(loglik, 1, max) 
	
	l.h <- exp(loglik - max.lik)

	result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
	return(result)
	
}

"likelihood.poisson" <- function(y, lambda) {

	N <- nrow(y)
	K <- ncol(lambda)
	nst <- nrow(lambda)
		
	y <- matrix(y, nrow = N, ncol = K)

	if(nst == 1) {
		lambda <- t(apply(matrix(1, nrow = N, ncol = K), 1, "*", lambda))
	}
	lambda <- apply(lambda, c(1,2), max, 0.0001, na.rm = TRUE)

	loglik <- y * log(lambda) - lambda - lgamma(y + 1)

	max.lik <- apply(loglik, 1, max, na.rm = TRUE)

	l.h <- exp(apply(loglik, 2, "-", max.lik)) 
	
	result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
	return(result) 
}

"likelihood.binomial" <- function(y, T, p) {

	N <- nrow(y)
	K <- ncol(p)
	nst <- nrow(T)
	
	y <- matrix(y, nrow = N, ncol = K)

	if(nst == 1) {
		T <- matrix(T, nrow = N, ncol = K)
	}	
	else {
		T <- t(apply(matrix(1, nrow = N, ncol = K), 1, "*", T))
	}
	
	loglik <- lgamma(T + 1) - lgamma(T - y + 1) - lgamma(y + 1)
	loglik <- loglik + t((apply(y, 1, "*", log(p)))) + t((apply(T - y, 1, "*", log(1 - p))))

	max.lik <- apply(loglik, 1, max, na.rm = TRUE)
	
	l.h <- exp(apply(loglik, 2, "-", max.lik))

	result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
	return(result)
}

"likelihood.normult" <- function(y, mu, sigmainv, logdet) {
	
	N <- nrow(y)
	r <- ncol(y)
	K <- dim(sigmainv)[3]
	loglik <- matrix(0, nrow = N, ncol = K)
	loglik1 <- -.5 * r * log(2*pi)

	for(k in 1:K) {
		eps <- t(apply(y, 1, "-",  t(matrix(mu[,k]))))
		loglik[, k] <- loglik1 + .5 * logdet[k] - .5 * apply(eps %*% sigmainv[,,1] * eps, 1, sum)		
	}

	maxlik <- t(apply(loglik, 1, max))
	l.h <- exp(apply(loglik, 2, "-", maxlik))

	results <- list(lh = l.h, maxl = matrix(maxlik), llh = loglik)
	return(results)
}

"likelihood.studmult" <- function(y, mu, sigmainv, logdet, df) {

	N <- nrow(y)
	K <- ncol(mu)
	r <- ncol(y)
 
	loglik <- matrix(0, nrow = N, ncol = K)
	
	for(k in 1:K) {
		eps <- t(apply(y, 1, "-", t(matrix(mu[,k]))))
		err <- apply(eps %*% sigmainv[,,k] * eps, 1, sum)
		loglik[, k] <- lgamma((df[k] + r)/2) - lgamma(df[k]/2) + .5 * logdet[k] - .5 * r * log(df[k] * pi) 
	        loglik[, k] <- loglik[, k] - ((df[k] + r)/2) * log(1 + err/df[k])
	}		

	maxlik <- t(apply(loglik, 1, max))
	l.h <- exp(apply(loglik, 2, "-", maxlik))

	result <- list(lh = l.h, maxl = matrix(maxlik), llh = loglik)
	return(result)
}
