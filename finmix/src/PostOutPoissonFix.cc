#include "PostOutPoissonFix.h"

PostOutPoissonFix::PostOutPoissonFix (const Rcpp::List& list) 
{
	Rcpp::List tmpList((SEXP) list["par"]);
	Rcpp::NumericMatrix tmpA((SEXP) tmpList["a"]);
	Rcpp::NumericMatrix tmpB((SEXP) tmpList["b"]);
	const unsigned int M = tmpA.nrow();
	const unsigned int K = tmpA.ncol();	
	a = new arma::mat(tmpA.begin(), M, K, false, true);
	b = new arma::mat(tmpB.begin(), M, K, false, true);
}

void PostOutPoissonFix::store (const unsigned int& m,
				const PriorPoissonFix& hyperPar) 
{
	(*a).row(m) = hyperPar.aPost;
	(*b).row(m) = hyperPar.bPost;
}
