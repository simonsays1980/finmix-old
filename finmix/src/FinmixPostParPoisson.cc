#include <RcppArmadillo.h>	
#include "FinmixPostParPoisson.h"
// post in R contains a list par
FinmixPostParPoisson::FinmixPostParPoisson (const Rcpp::List& postpar) 
{
	Rcpp::NumericMatrix tmpA(postpar["a"]);
	Rcpp::NumericMatrix tmpB(postpar["b"]);
	unsigned int M = tmpA.rows();
	unsigned int K = tmpA.cols();
	aM(tmpA.begin(), M, K, false, true);
	bM(tmpB.begin(), M, K, false, true);
}

