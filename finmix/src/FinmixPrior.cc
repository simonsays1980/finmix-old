/*
 * FinmixPrior.cc
 * 
 * Definition of FinmixPrior struct to contain 
 * finmix 'prior' S4-object 
 * Also a converter from Rcpp S4-object to
 * C++-struct is provided.
 *
 *  author: Lars Simon Zehnder
 * package: finmix (1.0.0)
 * created: 19 Feb. 2013
 */

#include <RcppArmadillo.h>
#include "FinmixPrior.h"

FinmixPrior::FinmixPrior(const Rcpp::S4& classS4) {

	/**
	 * As the parameters 'hyper_par' can have very differing implementations 
	 * they are kept in a general Rcpp::List object which is then 
	 * decomposed in main code.
	 */ 
	hyperParList = Rcpp::as<Rcpp::List>(classS4.slot("par"));
	weightV      = Rcpp::as<arma::rowvec>(classS4.slot("weight"));
	hier 	     = Rcpp::as<bool>(classS4.slot("hier"));
	priorType    = Rcpp::as<std::string>(classS4.slot("type"));
}

FinmixPrior::~FinmixPrior(){}
