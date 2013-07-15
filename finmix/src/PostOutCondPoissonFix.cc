#include "PostOutCondPoissonFix.h"

PostOutCondPoissonFix::PostOutCondPoissonFix (const Rcpp::List& list) 
{
	Rcpp::List tmpList((SEXP) list["par"]);
	Rcpp::NumericMatrix tmpA((SEXP) tmpList["a"]);
	Rcpp::NumericMatrix tmpB((SEXP) tmpList["b"]);
	Rcpp::NumericMatrix tmpCond((SEXP) tmpList["cond"]);
	const unsigned int M = tmpA.nrow();
	const unsigned int K = tmpA.ncol();	
	a = new arma::mat(tmpA.begin(), M, K, false, true);
	b = new arma::mat(tmpB.begin(), M, K, false, true);
	cond = new arma::mat(tmpCond.begin(), M, K, false, true);
}

void PostOutCondPoissonFix::store (const unsigned int& m,
	const PriorCondPoissonFix& hyperPar) 
{
	(*a).row(m) = hyperPar.aPost;
	(*b).row(m) = hyperPar.bPost;
	(*cond).row(m) = hyperPar.cond;
}
