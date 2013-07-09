#ifndef FINMIXPARPOISSON_H
#define FINMIXPARPOISSON_H

#include <RcppArmadillo.h>

class FinmixParPoisson {
	public: 
		arma::mat lambda;
		FinmixParPoisson(const Rcpp::List &par);
		~FinmixParPoisson();
};

#endif
