#ifndef PAROUTPOISSON_H
#define PAROUTPOISSON_H

#include <RcppArmadillo.h>
#include "ParPoissonFix.h"

class ParOutPoisson {
	public:
		arma::mat* lambda;

		ParOutPoisson () {}		
		ParOutPoisson (const Rcpp::List&); 
		void store(const unsigned int&, const ParPoissonFix&);
};
#endif
