#ifndef FINMIXPOSTPARPOISSON_H
#define FINMIXPOSTPARPOISSON_H

#include <RcppArmadillo.h>

class FinmixPostParPoisson {
	public:
		arma::mat aM;
		arma::mat bM;
		FinmixPostParPoisson (const Rcpp::List& postpar);
		~FinmixPostParPoisson();
};

#endif
