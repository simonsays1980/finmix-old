#include <RcppArmadillo.h>
#include "FinmixHierPoisson.h"

FinmixHierPoisson::FinmixHierPoisson (const Rcpp::List& hyper)
{
	Rcpp::NumericVector tmpB(hyper["b"]);
	unsigned int M = tmpB.nrows();
	hyper(tmpB.begin(), M, false, true);
}
