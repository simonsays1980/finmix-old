#ifndef PRIORPOISSONIND_H
#define PRIORPOISSONIND_H

#include <RcppArmadillo.h>
#include "PriorPoissonFix.h"

class ParPoissonInd;
class PriorPoissonInd : virtual public PriorPoissonFix {
	public:
		arma::rowvec weightStart;
		arma::rowvec weightPost;
		
		PriorPoissonInd (const FinmixPrior&);
		virtual ~PriorPoissonInd () {}
		virtual void update (const unsigned int&,
			const arma::mat&, arma::ivec&,
			const ParPoissonInd&);
};
#endif
