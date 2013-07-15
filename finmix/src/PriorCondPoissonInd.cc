#include "PriorCondPoissonInd.h"
#include "posterior.h"

PriorCondPoissonInd::PriorCondPoissonInd (const FinmixPrior& prior) :
	PriorPoissonFix(prior),
	PriorPoissonInd(prior),
	PriorCondPoissonFix(prior) {}

void PriorCondPoissonInd::update (const unsigned int& K, const arma::mat& y,
	arma::ivec& S, const ParPoissonInd& par)
{
	PriorPoissonInd::update(K, y, S, par);
}

void PriorCondPoissonInd::updateHier(const ParPoissonFix& par) 
{
	PriorCondPoissonFix::updateHier(par);
}
