#ifndef POSTOUTPOISSONFIX_H
#define POSTOUTPOISSONFIX_H

#include <RcppArmadillo.h>
#include "PriorPoissonFix.h"

class PostOutPoissonFix {
	public:
		arma::mat *a;
		arma::mat *b;

		PostOutPoissonFix () {}		
		PostOutPoissonFix (const Rcpp::List&);
		~PostOutPoissonFix () {}
		void store (const unsigned int&, 
			const PriorPoissonFix&);
};
#endif
