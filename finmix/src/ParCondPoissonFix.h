#ifndef PARCONDPOISSONFIX_H
#define PARCONDPOISSONFIX_H

#include <RcppArmadillo.h>
#include "ParPoissonFix.h"
#include "FinmixModel.h"
#include "PriorCondPoissonFix.h"
#include "distributions.h"

class ParCondPoissonFix : virtual public ParPoissonFix {
	public: 
		ParCondPoissonFix (const bool&, 
				const FinmixModel&);
		virtual ~ParCondPoissonFix () {}
		void update (PriorCondPoissonFix&);
};
#endif
