#include "PostOutCondPoissonInd.h"

PostOutCondPoissonInd::PostOutCondPoissonInd (const Rcpp::List& list) :
	PostOutCondPoissonFix(list) 
{	
	Rcpp::NumericMatrix tmpWeight((SEXP) list["weight"]);
	const unsigned int M = tmpWeight.nrow();
	const unsigned int K = tmpWeight.ncol();
	weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

void PostOutCondPoissonInd::store (const unsigned int& m,
	const PriorCondPoissonInd& hyperPar) 
{
	PostOutCondPoissonFix::store(m, hyperPar);
	(*weight).row(m) = hyperPar.weightPost;
}
