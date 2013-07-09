#ifndef PRIORPOISSONFIX_H
#define PRIORPOISSONFIX_H

#include <RcppArmadillo.h>
#include "FinmixPrior.h"

class ParPoissonFix;
class PriorPoissonFix {
	public:
		arma::rowvec aStart;
		arma::rowvec bStart;
		arma::rowvec aPost;
		arma::rowvec bPost;
		const bool HIER;
		double g;
		double G;

		PriorPoissonFix ();	
		PriorPoissonFix	(const FinmixPrior&);
		virtual	~PriorPoissonFix () {} 
		virtual void update (const unsigned int&, 
			const arma::mat&, arma::ivec&,
			const ParPoissonFix&);
		virtual void updateHier(const ParPoissonFix&);
};
#endif
