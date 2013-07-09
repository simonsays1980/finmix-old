/* 
 * Methods for calculating the log likelihood of used
 * priors. 
 * This class contrasts from 'likelihood.h' by not working 
 * with the 'liklist' struct.
 * 
 * author : Lars Simon Zehnder
 * package: finmix
 * created: 02/01/2013
 * changed: 25/05/2013
 */ 
#ifndef PRIORLIKELIHOOD_H
#define PRIORLIKELIHOOD_H

#include <RcppArmadillo.h>
#include <R.h> 			// to interface with R
#include <Rmath.h> 		// for using internal R C-functions
#include "likelihood.h"
/**
 * Evaluates the prior likelihood for the 
 * weights. This function is used in every
 * Gibbs sampling. 
 * 
 */
inline double 
priormixlik_dirichlet(const arma::rowvec &eta, 
			const arma::rowvec &prior_par) 
{
	unsigned int K = eta.n_elem;
	double priormixlik = 0.0;	
	/* evaluate dirichlet loglik */
	priormixlik = R::lgammafn(arma::accu(prior_par));
	for(unsigned int k = 0; k < K; ++k) {
		priormixlik += (prior_par(k) - 1) * std::log(eta(k));
			priormixlik -= R::lgammafn(prior_par(k));		 
	}
	return priormixlik;
}

inline double
priormixlik_poisson(const arma::rowvec& lambda, 
			const arma::rowvec& prior_parA,
			const arma::rowvec& prior_parB, 
			const bool &hier, const double &g, 
			const double &G) {

	unsigned int K = lambda.n_elem;
	double priormixlik = 0.0;
	if(!hier) {
		/*double scale = 0.0;
		for(unsigned int k = 0; k < K; ++k) {
			scale = 1.0/prior_parB(k);
			priormixlik += R::dgamma(lambda(k), prior_parA(k),
scale, 1);
		}*/
		priormixlik = likelihood_gamma(lambda, prior_parA(0), prior_parB(0));
	}
	else { /* hierarchical prior */
		double gN = g + K * prior_parA(0); // prior_parA must be the start value.
		double GN = G + arma::accu(lambda);
		double b = gN/GN;
		double scale = 1.0/b;
		/* step 1: log likelihood of prior */
	 	for(unsigned int k = 0; k < K; ++k) {
			priormixlik += R::dgamma(lambda(k), prior_parA(k), scale, 1);
		}
			/**
		 * step 2: log likelihood of hyperprior with start 
		 * values of hyper parameters 
		 */
		scale = 1.0/G;
		priormixlik += R::dgamma(b, g, scale, 1);
	
		/**
		 * step 3: log likelihood of hyperprior with updated 
		 * hyper parameters 
 		 */
		scale = 1.0/GN;
		priormixlik -= R::dgamma(b, gN, scale, 1);

	}
	return priormixlik;	
}
#endif
