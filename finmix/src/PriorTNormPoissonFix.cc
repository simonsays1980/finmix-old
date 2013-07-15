#include "PriorTNormPoissonFix.h"

PriorTNormPoissonFix::PriorTNormPoissonFix () {}


PriorTNormPoissonFix::PriorTNormPoissonFix (const FinmixPrior& prior) :
	sStart(Rcpp::as<arma::rowvec>((SEXP) prior.par["s"])) {}

void PriorTNormPoissonFix::update (const unsigned int& K, const arma::mat& y,
			arma::ivec& S, const ParPoissonFix& par)  {}

void PriorTNormPoissonFix::updateHier(const ParPoissonFix& par) {}

