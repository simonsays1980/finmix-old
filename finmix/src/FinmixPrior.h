/*
 * FinmixPrior.h
 * 
 * Definition of FinmixModel struct to contain 
 * finmix 'prior' S4-object 
 * Also a converter from Rcpp S4-object to
 * C++-struct is provided.
 *
 *  author: Lars Simon Zehnder
 * package: finmix (1.0.0)
 * created: 19 Feb. 2013
 */
#ifndef FINMIXPRIOR_H
#define FINMIXPRIOR_H

#include <RcppArmadillo.h>
#include <string>

class FinmixPrior {
	
	public:
		Rcpp::List par;
		arma::rowvec weight;
		
		std::string type;
		bool hier;
	
		/* ctor */ 
		FinmixPrior(const Rcpp::S4& classS4);
};
#endif
