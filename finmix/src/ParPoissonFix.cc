#include "ParPoissonFix.h"

ParPoissonFix::ParPoissonFix (const bool& STARTPAR, 
		const FinmixModel& model) : lambda(model.K) 
{
	if(STARTPAR && model.K > 1) {
		arma::rowvec tmp = Rcpp::as<arma::rowvec>
				((SEXP) model.par["lambda"]);
		lambda = tmp;
	}
} 

/*ParPoissonFix::ParPoissonFix(const FinmixModel& model) : 
	lambda(Rcpp::as<arma::rowvec>((SEXP) model.par["lambda"])){}
*/ 
void ParPoissonFix::update (const PriorPoissonFix& hyperPar) 
{
	lambda = rgammaprod(hyperPar.aPost, hyperPar.bPost);
}
