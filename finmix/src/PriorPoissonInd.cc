#include "PriorPoissonInd.h"
#include "ParPoissonInd.h"
#include "posterior.h"

PriorPoissonInd::PriorPoissonInd (const FinmixPrior& prior) :
	PriorPoissonFix(prior),
	weightStart(prior.weight),
	weightPost(prior.weight) {}

void PriorPoissonInd::update (const unsigned int& K, const arma::mat& y,
	arma::ivec& S, const ParPoissonInd& par) 

{
	PriorPoissonFix::update(K, y, S, par);
	weightPost = posterior_multinomial(K, S, weightStart);	
}
