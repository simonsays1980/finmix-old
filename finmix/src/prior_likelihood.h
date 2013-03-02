/* 
 * Methods for calculating the log likelihood of used
 * priors. 
 * This class contrasts from 'likelihod.h' by not working 
 * with the 'liklist' struct.
 * 
 * author : Lars Simon Zehnder
 * package: finmix
 * created: 02/01/2013
 */ 

#include <RcppArmadillo.h>
#include <Rmath.h>   // for use of R internal functions

double priormixlik_poisson(const arma::rowvec &lambda, 
			const arma::mat &prior_par,
			const bool &hier, const double &g, 
			const double &G) {

	unsigned int K = lambda.n_elem;
	double priormixlik = 0.0;
	if(!hier) {
		double scale = 0.0;
		for(unsigned int k = 0; k < K; ++k) {
			scale = 1.0/prior_par(k, 1);
			priormixlik += R::dgamma(lambda(k), prior_par(k, 0), scale, 1);
		}
	}
	else { /* hierarchical prior */
		double gN = g + K * prior_par(0, 0); // prior_par must be the start value.
		double GN = G + arma::sum(lambda);
		double b = gN/GN;
		double scale = 1.0/b;
		/* step 1: log likelihood of prior */
	 	for(unsigned int k = 0; k < K; ++k) {
			priormixlik += R::dgamma(lambda(k), prior_par(k, 0), scale, 1);
		}	
		
		/**
		 * step 2: log likelihood of hyperprior with start 
		 * values of hyperparameters 
		 */
		priormixlik += K * R::dgamma(b, g, G, 1);
	
		
		/**
		 * step 3: log likelihood of hyperprior with updated 
		 * hyper parameters 
 		 */
		priormixlik -= K * R::dgamma(b, gN, GN, 1);
	}
	
	return priormixlik;
}
