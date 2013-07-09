#include "LogPoissonInd.h"

LogPoissonInd::LogPoissonInd () : LogPoissonFix(),
	cdpost(0.0), entropy(0.0), maxcdpost(0.0) {}

void LogPoissonInd::update (const unsigned int& K, 
	const arma::mat& y, arma::ivec &S, const arma::mat& expos,
	const ParPoissonInd& par, const PriorPoissonInd& hyperPar)
{
	
	arma::mat lambdaM = arma::kron(expos, par.lambda);
	liklist lik = likelihood_poisson(y, lambdaM);
	DataClass dataC = classification(S, lik, par.weight);
	S = dataC.newS;
	mixlik = dataC.mixLik;
		/* Compute likelihood of mixture prior */
	mixprior = priormixlik_poisson(par.lambda,
		hyperPar.aStart, hyperPar.bStart,
		hyperPar.HIER, hyperPar.g, hyperPar.G);
	if(K > 1) {
		/* Compute likelihood of Dirichlet prior */
		mixprior += priormixlik_dirichlet(par.weight, 
			hyperPar.weightStart);
		cdpost = mixlik + mixprior + dataC.postS;
		entropy = dataC.entropy;
	}
	
}
