#include "ParCondPoissonFix.h"

ParCondPoissonFix::ParCondPoissonFix (const bool& STARTPAR, 
	const FinmixModel& model) : 
		ParPoissonFix::ParPoissonFix(true, model) {}

void ParCondPoissonFix::update (PriorCondPoissonFix& hyperPar) 
{
	const unsigned int K = lambda.n_elem;
	double tmp = 0.0;
	for(unsigned int k = 0; k < K; ++k) {
		tmp = arma::as_scalar(hyperPar.coef.row(k) 
			* lambda.t());
		tmp -= lambda(k);
		hyperPar.cond(k) = tmp;	
		lambda(k) = rggamma(hyperPar.aPost(k), 
			hyperPar.bPost(k), tmp);
	}
}
