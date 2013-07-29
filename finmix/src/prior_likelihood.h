/******************************************************************************
 *
 * Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is part of the R package finmix.
 *
 * finmix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundatio, either version 3 of the License, or
 * any later version.
 *
 * finmix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with finmix. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
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

inline double 
priormixlik_condpoisson (const arma::rowvec& lambda,
	const arma::rowvec& prior_parA, const arma::rowvec& prior_parB,
	const arma::rowvec& prior_cond, const bool& HIER, const double &g,
	const double &G)
{
	unsigned int K = lambda.n_elem;
	double priormixlik = 0.0;
	if(!HIER) {
		priormixlik = likelihood_ggamma(lambda, prior_parA, 
			prior_parB(0), prior_cond);
	} else {
		double gN = g + K * prior_parA(0); 	// prior_parA must be the start value
		double GN = G + arma::accu(lambda);
		double b = gN/GN;
		double scale = 1.0/b;
		/* step 1: log likelihood of prior */
		for(unsigned int k = 0; k < K; ++k) {
			priormixlik += R::dgamma(lambda(k) - prior_cond(k), 
				prior_parA(k), scale, 1);
		}
		/**
		 * step 2: log likelihood of hyperprior with 
		 * 	start values of hyper parameters.
		 */
		scale = 1.0/G;
		priormixlik += R::dgamma(b, g, scale, 1);
		/**
		 * step 3: log likelihood of hyperprior with 
		 * 	updated hyper parameters.
		 */
		scale = 1.0/GN;
		priormixlik -= R::dgamma(b, gN, scale, 1);
	}
	return priormixlik;
} 
#endif
