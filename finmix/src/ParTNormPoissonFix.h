#ifndef PARTNORMPOISSONFIX_H
#define PARTNORMPOISSONFIX_H

#include <RcppArmadillo.h>
#include "ParPoissonFix.h"
#include "FinmixModel.h"
#include "PriorTNormPoissonFix.h"
#include "distributions.h"

class ParTNormPoissonFix : virtual public ParPoissonFix {
	public: 
		ParTNormPoissonFix (const bool&, 
				const FinmixModel&);
		virtual ~ParTNormPoissonFix () {}
		void update (PriorTNormPoissonFix&);
};
#endif
