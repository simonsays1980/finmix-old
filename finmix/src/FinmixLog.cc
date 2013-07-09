#include <RcppArmadillo.h>
#include "FinmixLog.h"

FinmixLog::FinmixLog (const Rcpp::List &log)
{
	Rcpp::NumericVector tmpML(log["mixlik"]);
	Rcpp::NumericVector tmpMP(log["mixprior"]);
	Rcpp::NumericVector tmpCD(log["cdpost"]);
	unsigned int M = tmpML.nrow();
	mixlik(tmpML.begin(), M, false, true);
	mixprior(tmpMP.begin(), M, false, true);
	cdpost(tmpCD.begin(), M, false, true);
}
