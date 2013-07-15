#include "PriorCondPoissonFix.h"
#include "ParCondPoissonFix.h"

PriorCondPoissonFix::PriorCondPoissonFix () : PriorPoissonFix() {}


PriorCondPoissonFix::PriorCondPoissonFix (const FinmixPrior& prior) :
	PriorPoissonFix(prior),
	coef(Rcpp::as<arma::mat>((SEXP) prior.par["coef"])),
	cond(Rcpp::as<arma::rowvec>((SEXP) prior.par["a"])) {}

void PriorCondPoissonFix::update (const unsigned int& K, const arma::mat& y,
			arma::ivec& S, const ParPoissonFix& par)  
{
	PriorPoissonFix::update(K, y, S, par);
}

void PriorCondPoissonFix::updateHier(const ParPoissonFix& par)
{
	if (HIER) {
		const unsigned int K = par.lambda.n_elem;
		GetRNGstate();
		double gN = g + arma::sum(aStart);
		double GN = 0.0;
		double b = 0.0;
		for(unsigned int k = 0; k < K; ++k) {
			GN = G * std::pow(2, arma::sum(coef.row(k)) - 1)
				+ arma::sum(par.lambda);
			b = R::rgamma(gN, 1/GN);
			bStart(k) = b;
		}
		PutRNGstate();
	}
}

