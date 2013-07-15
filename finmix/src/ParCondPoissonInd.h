#ifndef PARCONDPOISSONIND_H
#define PARCONDPOISSONIND_H

#include <RcppArmadillo.h>
#include "ParPoissonInd.h"
#include "ParCondPoissonFix.h"
#include "PriorCondPoissonInd.h"

class ParCondPoissonInd : public ParPoissonInd,
	public ParCondPoissonFix {
	public:	
		ParCondPoissonInd (const bool&, 
			const FinmixModel&);
		virtual ~ParCondPoissonInd () {}
		void update (PriorCondPoissonInd&);
};
#endif
