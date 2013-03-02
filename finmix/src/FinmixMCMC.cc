/*
 * FinmixMCMC.cc
 * 
 * Definition of FinmixMCMC struct to contain 
 * finmix 'model' S4-object 
 * Also a converter from Rcpp S4-object to
 * C++-struct is provided.
 *
 *  author: Lars Simon Zehnder
 * package: finmix (1.0.0)
 * created: 19 Feb. 2013
 */

#include <RcppArmadillo.h>
#include "FinmixMCMC.h"

FinmixMCMC::FinmixMCMC(const Rcpp::S4& classS4) {

	burnIn 	  = Rcpp::as<unsigned int>(classS4.slot("burnin"));
	m  	  = Rcpp::as<unsigned int>(classS4.slot("M"));
	startPar  = Rcpp::as<bool>(classS4.slot("startpar"));
	storeS 	  = Rcpp::as<unsigned int>(classS4.slot("storeS"));
	storePost = Rcpp::as<bool>(classS4.slot("storepost"));
	ranPerm   = Rcpp::as<bool>(classS4.slot("ranperm"));	
}

FinmixMCMC::~FinmixMCMC(){}


