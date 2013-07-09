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
#ifndef FINMIXMODEL_H
#define FINMIXMODEL_H

#include <RcppArmadillo.h>
#include <string>

class FinmixModel {
	
	public:
		Rcpp::List par;
		arma::rowvec weight;
		arma::ivec T;
	
		bool indicFix; 
		unsigned int K;
		unsigned int r;
	
		/* ctor */ 
		FinmixModel(const Rcpp::S4& classS4);	
		/* dtor */
		~FinmixModel();
};
#endif
