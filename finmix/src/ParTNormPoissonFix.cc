#include "ParTNormPoissonFix.h"

ParTNormPoissonFix::ParTNormPoissonFix (const bool& STARTPAR, 
	const FinmixModel& model) : 
		ParPoissonFix::ParPoissonFix(true, model) {}

void ParTNormPoissonFix::update (PriorTNormPoissonFix& hyperPar) 
{
	const unsigned int K = lambda.n_elem;
	double tmp = 0.0;
	
	for(unsigned int k = 1; k < K - 1; ++k) {
		tmp = std::pow(lambda(k - 1) - lambda(k), 2)
			+ std::pow(lambda(k) - lambda(k + 1), 2);
	}
		
}
