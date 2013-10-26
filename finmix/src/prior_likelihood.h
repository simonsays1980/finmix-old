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

// =============================================================
// Dirichlet prior
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Computes the mixture log-likelihood for a 
 *          Dirichlet distribution. 
 * @par eta         weight parameters
 * @par prior_par   hyper parameters of the Dirichlet 
 *                  distribution
 * @return  Prior mixture log-likelihood
 * @detail  The log-lieklihood is computed via the density
 *          of the Dirichlet distribution evaluated at the 
 *          weights with hyper parameters determined in 
 *          prior_par
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
inline double 
priormixlik_dirichlet(const arma::rowvec &eta, 
			const arma::rowvec &prior_par) 
{
	unsigned int K      = eta.n_elem;
	double priormixlik  = 0.0;	
	/* Evaluate Dirichlet loglik */
	priormixlik = R::lgammafn(arma::accu(prior_par));
	for(unsigned int k = 0; k < K; ++k) {
		priormixlik += (prior_par(k) - 1) * std::log(eta(k));
        priormixlik -= R::lgammafn(prior_par(k));		 
	}
	return priormixlik;
}
// =============================================================
// Poisson prior
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Computes the prior mixture log-likelihood for a 
 *          Poisson distribution. 
 * @par lambda      Poisson parameters, 1 x K 
 * @par prior_parA  Gamma hyper shape parameters, 1 x K
 * @par prior_parB  Gamma hyper rate parameters, 1 x K
 * @par hier        boolean value indicating if a hierarchical
 *                  prior is used
 * @par g           Gamma hyper shape parameter of the hierar-
 *                  chical prior
 * @par G           Gamma hyper rate parameter of the hierar-
 *                  chical prior
 * @return  the prior mixture log-likelihood value
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
inline double
priormixlik_poisson(const arma::rowvec& lambda, 
			const arma::rowvec& prior_parA,
			const arma::rowvec& prior_parB, 
			const bool &hier, const double &g, 
			const double &G) {

	unsigned int K = lambda.n_elem;
	double priormixlik = 0.0;
	if(!hier) {
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

// =============================================================
// Binomial distribution
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Computes the prior mixture log-lieklihood for a 
 *          Binomial distribution. 
 * @par lambda      Poisson parameters, 1 x K 
 * @par prior_parA  Gamma hyper shape parameters, 1 x K
 * @par prior_parB  Gamma hyper rate parameters, 1 x K
 * @par hier        boolean value indicating if a hierarchical
 *                  prior is used
 * @par g           Gamma hyper shape parameter of the hierar-
 *                  chical prior
 * @par G           Gamma hyper rate parameter of the hierar-
 *                  chical prior
 * @return  the prior mixture log-likelihood value
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
inline
double priormixlik_binomial (const arma::rowvec& p,
        const arma::rowvec& prior_parA, 
        const arma::rowvec& prior_parB)
{
     const unsigned int K = p.n_elem;
     double priormixlik = 0.0;
     for (unsigned int k = 0; k < K; ++k) {
         priormixlik += (prior_parA(k) - 1.0) * std::log(p(k)) 
             + (prior_parB(k) - 1.0) * std::log(p(k));
         priormixlik -= R::lbeta(prior_parA(k), prior_parB(k));
     }
     return priormixlik;
}
#endif
