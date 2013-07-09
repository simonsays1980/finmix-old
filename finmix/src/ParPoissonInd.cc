#include "ParPoissonInd.h"

ParPoissonInd::ParPoissonInd (const bool& STARTPAR, 
	const FinmixModel& model) :
		ParPoissonFix(STARTPAR, model), 
		weight(model.K) 
{
	if(STARTPAR && model.K > 1) {
		weight = model.weight;
	}
}

void ParPoissonInd::update (const PriorPoissonInd& hyperPar)
{
	ParPoissonFix::update(hyperPar);
	weight = rdirichlet(hyperPar.weightPost);
}
