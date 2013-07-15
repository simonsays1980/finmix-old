#ifndef PARPOISSONFIX_H
#define PARPOISSONFIX_H

#include <RcppArmadillo.h>
#include "FinmixModel.h"
#include "PriorPoissonFix.h"
#include "distributions.h"

class ParPoissonFix {
	public: 
		arma::rowvec lambda;
		
		ParPoissonFix (const bool&, 
				const FinmixModel&);
		virtual ~ParPoissonFix () {}
		void update (const PriorPoissonFix&);
};
#endif
