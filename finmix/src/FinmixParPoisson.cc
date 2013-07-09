#include <RcppArmadillo.h>
#include "FinmixParPoisson.h"

FinmixParPoisson::FinmixParPoisson (const Rcpp::List& par) 
{
	Rcpp::NumericMatrix tmpPar(par["lambda"]);
	unsigned int M = tmpPar.rows();
	lambda(tmpPar.begin(), M, false, true);
}
