#ifndef PRIORCONDPOISSONIND_H
#define PRIORCONDPOISSONIND_H

#include <RcppArmadillo.h>
#include "PriorCondPoissonFix.h"
#include "PriorPoissonInd.h"

class PriorCondPoissonInd : public PriorCondPoissonFix,
	public PriorPoissonInd {
	public:	
		PriorCondPoissonInd (const FinmixPrior&);
		virtual ~PriorCondPoissonInd () {}
		virtual void update (const unsigned int&, 
			const arma::mat&, arma::ivec&,
			const ParPoissonInd&);
		virtual void updateHier (const ParPoissonFix&);
};
#endif
