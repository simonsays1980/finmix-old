#ifndef LOGCONDPOISSONFIX_H
#define LOGCONDPOISSONFIX_H

#include <RcppArmadillo.h>
#include "ParCondPoissonFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogCondPoissonFix {
	public:
		double mixlik;
		double mixprior;
		
		LogCondPoissonFix ();
		virtual ~LogCondPoissonFix () {}
		void update (const unsigned int&, const arma::mat&, 
				const arma::ivec&, const arma::mat& expos,
				const ParCondPoissonFix&, const PriorCondPoissonFix&);
};
#endif
