#ifndef FINMIXHIERPOISSON_H
#define FINMIXHIERPOISSON_H

#include <RcppArmadillo.h>

class FinmixHierPoisson {
	public:
		arma::vec hyper;
		FinmixHierPoisson(const Rcpp::List& hyper);
		~FinmixHierPoisson;
};

#endif
