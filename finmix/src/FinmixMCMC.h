/*
 * FinmixMCMC.h
 * 
 * Definition of FinmixMCMC struct to contain 
 * finmix 'mcmc' S4-object 
 * Also a converter from Rcpp S4-object to
 * C++-struct is provided.
 *
 *  author: Lars Simon Zehnder
 * package: finmix (1.0.0)
 * created: 19 Feb. 2013
 */
#ifndef FINMIXMCMC_H
#define FINMIXMCMC_H

#include <RcppArmadillo.h>

class FinmixMCMC {
	
	public:
		unsigned int burnIn;
		unsigned int M;
		bool startPar;
		unsigned int storeS;
		bool storePost;
		bool ranPerm;
	
		/* ctor */ 
		FinmixMCMC(const Rcpp::S4& classS4);	
		/* dtor */
		~FinmixMCMC();
};
#endif