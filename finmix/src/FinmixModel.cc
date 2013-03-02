/*
 * FinmixModel.cc
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
#include "FinmixModel.h"

FinmixModel::FinmixModel(const Rcpp::S4& classS4) {

	/**
	 * As the parameters 'par' can have very differing implementations 
	 * they are kept in a general Rcpp::List object which is then 
	 * decomposed in main code.
	 *
	 */ 
	parList  = Rcpp::as<Rcpp::List>(classS4.slot("par"));
	weightV  = Rcpp::as<arma::rowvec>(classS4.slot("weight"));

	/** 
	 * 'T' can be 'NA' if the model has no repetitions (Binomial)
	 * in this case this attribute is set to an 1 x 1 vector
	 *
	 */
	tV	 = Rcpp::as<arma::uvec>(classS4.slot("T"));
	indicFix = Rcpp::as<bool>(classS4.slot("indicfix"));
	k 	 = Rcpp::as<unsigned int>(classS4.slot("K"));
	r 	 = Rcpp::as<unsigned int>(classS4.slot("r"));	
}

FinmixModel::~FinmixModel(){}
