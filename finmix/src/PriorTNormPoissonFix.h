#ifndef PRIORTNORMPOISSONFIX_H
#define PRIORTNORMPOISSONFIX_H

#include <RcppArmadillo.h>
#include "FinmixPrior.h"

class PriorTNormPoissonFix {
	public:
		arma::rowvec sStart;

		PriorTNormPoissonFix ();	
		PriorTNormPoissonFix	(const FinmixPrior&);
		virtual	~PriorTNormPoissonFix () {} 
		virtual void update (const unsigned int&, 
			const arma::mat&, arma::ivec&,
			const ParPoissonFix&);
		virtual void updateHier(const ParPoissonFix&);
};
#endif
