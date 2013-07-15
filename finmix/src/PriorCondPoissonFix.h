#ifndef PRIORCONDPOISSONFIX_H
#define PRIORCONDPOISSONFIX_H

#include <RcppArmadillo.h>
#include "PriorPoissonFix.h"
#include "FinmixPrior.h"

class PriorCondPoissonFix : virtual public PriorPoissonFix {
	public:
		arma::mat coef;
		arma::rowvec cond;

		PriorCondPoissonFix ();	
		PriorCondPoissonFix	(const FinmixPrior&);
		virtual	~PriorCondPoissonFix () {} 
		virtual void update (const unsigned int&, 
			const arma::mat&, arma::ivec&,
			const ParPoissonFix&);
		virtual void updateHier(const ParPoissonFix&);
};
#endif
