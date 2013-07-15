#ifndef PARPOISSONIND_H
#define PARPOISSONIND_H

#include <RcppArmadillo.h>
#include "ParPoissonFix.h"
#include "PriorPoissonInd.h"

class ParPoissonInd : virtual public ParPoissonFix {
	public:
		arma::rowvec weight;
		
		ParPoissonInd (const bool&, 
			const FinmixModel&);
		virtual ~ParPoissonInd () {}
		void update (const PriorPoissonInd&);
};
#endif
