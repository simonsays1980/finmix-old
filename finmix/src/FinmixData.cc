/**
* 'exp' and 'T' can be 'NA' if the model
* has no exposures (Poisson) or repetitions (Binomial)
* in this case these attributes are set to an 1 x 1 vector
*
*/

#include "FinmixData.h"

FinmixData::FinmixData (const Rcpp::S4& classS4) :
	y(Rcpp::as<arma::mat>((SEXP) classS4.slot("y"))),
	S(Rcpp::as<arma::ivec>((SEXP) classS4.slot("S"))),
	expos(Rcpp::as<arma::vec>((SEXP) classS4.slot("exp"))),
	T(Rcpp::as<arma::ivec>((SEXP) classS4.slot("T"))),
	dataType(Rcpp::as<std::string>((SEXP) classS4.slot("type"))),
	N(Rcpp::as<unsigned int>((SEXP) classS4.slot("N"))),
	r(Rcpp::as<unsigned int>((SEXP) classS4.slot("r"))) {}		
	
