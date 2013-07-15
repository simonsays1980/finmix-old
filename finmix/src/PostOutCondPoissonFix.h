#ifndef POSTOUTCONDPOISSONFIX_H
#define POSTOUTCONDPOISSONFIX_H

#include <RcppArmadillo.h>
#include "PriorCondPoissonFix.h"

class PostOutCondPoissonFix {
	public:
		arma::mat* a;
		arma::mat* b;
		arma::mat* cond;
	
		PostOutCondPoissonFix () {}		
		PostOutCondPoissonFix (const Rcpp::List&);
		~PostOutCondPoissonFix () {}
		void store (const unsigned int&, 
			const PriorCondPoissonFix&);
};
#endif
