#include "PostOutPoissonInd.h"

PostOutPoissonInd::PostOutPoissonInd (const Rcpp::List& list) :
	PostOutPoissonFix(list) 
{	
	Rcpp::NumericMatrix tmpWeight((SEXP) list["weight"]);
	const unsigned int M = tmpWeight.nrow();
	const unsigned int K = tmpWeight.ncol();
	weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

void PostOutPoissonInd::store (const unsigned int& m,
	const PriorPoissonInd& hyperPar) 
{
	PostOutPoissonFix::store(m, hyperPar);
	(*weight).row(m) = hyperPar.weightPost;
}
