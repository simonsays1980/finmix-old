#include "LogCondPoissonFix.h"

LogCondPoissonFix::LogCondPoissonFix () : mixlik(0.0), 
				mixprior(0.0) {}

void LogCondPoissonFix::update (const unsigned int& K, const arma::mat& y, 
			const arma::ivec& S, const arma::mat& expos, 
			const ParCondPoissonFix& par, 
			const PriorCondPoissonFix& hyperPar) 
{
	arma::mat lambdaM = arma::kron(expos, par.lambda);
	liklist lik = likelihood_poisson(y, lambdaM);
	DataClass dataC = classification_fix(K, S, lik);
	mixlik = arma::sum(dataC.logLikCd);
	/* Compute likelihood of mixture prior */
	mixprior = priormixlik_condpoisson(par.lambda, hyperPar.aStart, 
			hyperPar.bStart, hyperPar.cond,	hyperPar.HIER, 
			hyperPar.g, hyperPar.G);
}
