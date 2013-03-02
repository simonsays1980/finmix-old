/*
 * FinmixModel.h
 * 
 * Definition of FinmixModel struct to contain 
 * finmix 'model' S4-object 
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

class FinmixModel {
	
	public:
		Rcpp::List parList;
		arma::rowvec weightV;
		arma::uvec tV;
	
		std::string dataType;
		bool indicFix; 
		unsigned int k;
		unsigned int r;
	
		/* constructor */ 
		FinmixModel(const Rcpp::S4& classS4);	
		/* dtor */
		~FinmixModel();
};

