#include "LogPoissonFix.h"

LogPoissonFix::LogPoissonFix () : mixlik(0.0), 
				mixprior(0.0) {}

void LogPoissonFix::update (const unsigned int& K, const arma::mat& y, 
			const arma::ivec& S, const arma::mat& expos, 
			const ParPoissonFix& par, 
			const PriorPoissonFix& hyperPar) 
{
	arma::mat lambdaM = arma::kron(expos, par.lambda);
	liklist lik = likelihood_poisson(y, lambdaM);
	DataClass dataC = classification_fix(K, S, lik);
	mixlik = arma::sum(dataC.logLikCd);
	/* Compute likelihood of mixture prior */
	mixprior = priormixlik_poisson(par.lambda, hyperPar.aStart, hyperPar.bStart,
			hyperPar.HIER, hyperPar.g, hyperPar.G);
}
