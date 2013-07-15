#include "ParCondPoissonInd.h"

ParCondPoissonInd::ParCondPoissonInd (const bool& STARTPAR, 
	const FinmixModel& model) :
		ParPoissonFix(STARTPAR, model),
		ParPoissonInd(STARTPAR, model),
		ParCondPoissonFix(STARTPAR, model) {}

void ParCondPoissonInd::update (PriorCondPoissonInd& hyperPar)
{
	ParCondPoissonFix::update(hyperPar);
	weight = rdirichlet(hyperPar.weightPost);
}
