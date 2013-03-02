/*
 * MCMC sampler for poisson distributions
 * author: Lars Simon Zehnder
 * package: finmix 
 * created: 01/31/2013
 *
 */

#ifndef _finmix_MCMC_POISSON_CC
#define _finmix_MCMC_POISSON_CC

#include <RcppArmadillo.h>		// C++ linear algebra library
// [[Rcpp::depends(RcppArmadillo)]]

//#include "parprior_poisson.h"		//
#include "FinmixData.h"			// FinmixData class
#include "FinmixModel.h"		// FinmixModel class
#include "FinmixPrior.h"		// FinmixPrior class
#include "FinmixMCMC.h"			// FinmixMCMC class
#include "distributions.h"              // several distributions not in R
#include "posterior.h"			// calculate posterior update a_k, b_k, e_k of poisson
#include "DataClass.h"         		// for dataclass structs 
#include "prior_likelihood.h"  		// for log likelihoods of priors

#include <R.h>             		// needed to use Rprintf()
#include <R_ext/Utils.h>   		// needed to allow user interrupts 
#include <Rdefines.h>      		// for sharing code with S at any stage
#include <Rinternals.h>   		// defines handling of R objects in C


void mcmc_poisson_impl (const FinmixData<arma::vec>& finData, const FinmixModel& finModel, 
			const FinmixPrior& finPrior, const FinmixMCMC& finMCMC) {
		//TODO: Check mcmc paramaters from 'mcmc' object!!! 
		//TODO: Check for indicfix!!!
		//TODO: Chec for starting with parameters and starting with allocations
		const unsigned int K = finModel.k;						// number of components
		const bool HIER      = finPrior.hier;						// hierarchical prior
		const bool INDICFIX  = finModel.indicFix;
		
		const arma::rowvec hyperWeightStart(finPrior.weightV);				// starting hyper-
												// parameters for 
												// for the weights
		arma::vec  yM	     = finData.yM; 
		arma::uvec sV 	     = finData.sV;
		/**
		 * get hyper parameters from FinmixPrior object 
		 * 'parList' contains different number of elements:
		 * hier: TRUE
		 * 	a : an K x 1 matrix 
		 * 	b : an K x 1 matrix
		 * 	g : numeric 
		 * 	G : numeric
		 * 
		 * hier FALSE
		 * 	a : an K x 1 matrix 
		 * 	b : an K x 1 matrix
		 * 
		 */
		arma::mat hyperParStart(K, 2);

		hyperParStart.col(0) = Rcpp::as<arma::vec>(finPrior.hyperParList["a"]); 	// hyper parameter 'a' for
												// lambda
		hyperParStart.col(1) = Rcpp::as<arma::vec>(finPrior.hyperParList["b"]);		// hyper parameter 'b' for
												// lambda
		arma::rowvec eta(K); 					 	                // actual weights in 
												// step 'i', 1 x K vector
		arma::rowvec mixPar(K);		 						// actual component 
												// parameters in step 'i'
												// 1 x K vector 
		arma::mat hyperPar(K, 2); 						      	// matrix containing:
						 						// 'a' in first column
						 						// 'b' in second column
												// K x 2 matrix    
		arma::rowvec hyperWeight(K);       						// vector for the actual 
												// hyperparameters
											 	// for the weights in 
												// step 'i'
		double mixLik 	   = 0.0; 							 	// mixture likelihood 
		double priorMixLik = 0.0; 						 	// likelihood of the prior
		double g 	   = 0.0;									// shape of Gamma
		double G 	   = 0.0;									// rate of Gamma
		if(HIER) {
			g = finPrior.hyperParList("g");
			G = finPrior.hyperParList("G");
		}
		double gN 	   = 0.0;			 					// containing the updated
												// shape of hier Gamma
		double GN 	   = 0.0;			 					// containing the updated
												// rate of hier Gamma
		double b 	   = 0.0; 		 						// containing the updated
												// hyperparameter b 
												// (the rate) for lambda

		/* update parameters */ 

			/* weights parameters */
				/**
			 	 * update dirichlet parameters
				 * use the scheme:
				 * e_k(S) = e_0 + N_k(S)
				 */
		hyperWeight 	   = posterior_multinomial(K, sV, hyperWeightStart);
		hyperWeight.print("hyperWeight:");
				/* update weights of mixture */
		eta 		   = rdirichlet(hyperWeight);
		eta.print("eta:");
			/* component parameters */		
				/** 
				 * update parameters of prior for component parametes
				 * use with the following scheme:
				 * a_k(S) = a_0 + N_k(S) * mean(y)_k(S) 
				 * b_k(S) = b_0 + N_k(S)
				 */
		hyperParStart.print("hyperParStart:");
		hyperPar 	   = posterior_poisson(K, yM, sV, hyperParStart);
		hyperPar.print("hyperPar:");
				/**
				 * update component parameters
				 * the posterior is p(lambda_k|S, y) ~ G(a_k(S), b_k(S)) 
				 */
					//TODO: Check for input matrix!
		mixPar 		   = rgammaprod(hyperPar.col(0), hyperPar.col(1));
		mixPar.print("mixPar");
		Rprintf("HIER:%i\n", HIER);
		if(HIER) { // hierarchical prior for 'b'
			GetRNGstate(); // get RNG state from R
			// sample from G(g_0 + Ka_0, G_0 + sum lambda_k)
			gN = g + sum(hyperParStart.col(0));
			GN = G + sum(mixPar);

			b  = R::rgamma(gN, GN);

			PutRNGstate(); // close rng from R
			/* update random hyper parameter */
			hyperParStart.col(1).fill(b);			
		}

		/* update allocations and calculate mixture likelihood */
		if(!INDICFIX) {
 			DataClass dataC = dataclass_poisson(INDICFIX, yM, sV, eta, mixPar);
			sV		= dataC.newS;
			mixLik 		= dataC.mixLik;
		}
		else { /* fixed indicators */
			DataClass dataC = dataclass_poisson(INDICFIX, yM, sV, eta, mixPar);
			mixLik 		= arma::sum(dataC.logLikCd);			
		}
		
		/* compute likelihood of prior (gamma) */
		priorMixLik = priormixlik_poisson(mixPar, hyperParStart, HIER, g, G);
		Rprintf("priorMixLik: %15.9f\n");
}

RcppExport SEXP mcmc_poisson_cc(SEXP data_S4, SEXP model_S4, 
				SEXP prior_S4, SEXP mcmc_S4) {

	/* Convert S4-classes to C++-structs */
	Rcpp::S4 dataS4O(data_S4);
	Rcpp::S4 modelS4O(model_S4);
	Rcpp::S4 priorS4O(prior_S4);
	Rcpp::S4 mcmcS4O(mcmc_S4);
	FinmixData<arma::vec> finData = FinmixData<arma::vec>(dataS4O);
	FinmixModel finModel = FinmixModel(modelS4O);
	FinmixPrior finPrior = FinmixPrior(priorS4O);
	FinmixMCMC finMCMC = FinmixMCMC(mcmcS4O);
	mcmc_poisson_impl(finData, finModel, finPrior, finMCMC);	

	return Rcpp::wrap(1);	
}
#endif
