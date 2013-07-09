#ifndef LOGPOISSONFIX_H
#define LOGPOISSONFIX_H

#include <RcppArmadillo.h>
#include "ParPoissonFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogPoissonFix {
	public:
		double mixlik;
		double mixprior;
		
		LogPoissonFix ();
		virtual ~LogPoissonFix () {}
		void update (const unsigned int&, const arma::mat&, 
				const arma::ivec&, const arma::mat& expos,
				const ParPoissonFix&, const PriorPoissonFix&);
};
#endif
