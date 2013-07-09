#ifndef HIEROUTPOISSON_H
#define HIEROUTPOISSON_H

#include <RcppArmadillo.h>

class HierOutPoisson {
	public:
		arma::vec* b;

		HierOutPoisson () {}			
		HierOutPoisson (const Rcpp::List&);
		template <typename PriorParType>
		void store (const unsigned int& m, 
				const PriorParType&);
}; 

HierOutPoisson::HierOutPoisson (const Rcpp::List& list) 
{
	Rcpp::NumericVector tmpB((SEXP) list["b"]);
	const unsigned int M = tmpB.size();
	b = new arma::vec(tmpB.begin(), M, false, true);
}

template <typename PriorParType>
void HierOutPoisson::store (const unsigned int& m,
				const PriorParType& hyperPar) 
{
	(*b)(m) = hyperPar.bStart(0);
}
#endif
