#include "ParOutPoisson.h"

ParOutPoisson::ParOutPoisson(const Rcpp::List& list) 
{
	Rcpp::NumericMatrix tmpLambda((SEXP) list["lambda"]);
	const unsigned int M = tmpLambda.nrow();
	const unsigned int K = tmpLambda.ncol();
	lambda = new arma::mat(tmpLambda.begin(), M, K, false, true); 
} 

void ParOutPoisson::store (const unsigned int& m, const ParPoissonFix& par)
{
	(*lambda).row(m) = par.lambda;
}
