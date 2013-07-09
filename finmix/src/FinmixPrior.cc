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

#include "FinmixPrior.h"

FinmixPrior::FinmixPrior(const Rcpp::S4& classS4):
	par(Rcpp::as<Rcpp::List>((SEXP) classS4.slot("par"))),
	weight(Rcpp::as<arma::rowvec>((SEXP) classS4.slot("weight"))),
	type(Rcpp::as<std::string>((SEXP) classS4.slot("type"))),
	hier(Rcpp::as<bool>((SEXP) classS4.slot("hier"))) {}
