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

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <string>

class FinmixPrior {
	
	public:
		Rcpp::List hyperParList;
		arma::rowvec weightV;
		
		std::string priorType;
		bool hier;
	
		/* ctor */ 
		FinmixPrior(const Rcpp::S4& classS4);	
		/* dtor */
		~FinmixPrior();
};

