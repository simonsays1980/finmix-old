#ifndef LOGPOISSONIND_H
#define LOGPOISSONIND_H

#include <RcppArmadillo.h>
#include "LogPoissonFix.h"
#include "ParPoissonInd.h"
#include "PriorPoissonInd.h"

class LogPoissonInd : public LogPoissonFix {
	public:
		double cdpost;
		double entropy;
		double maxcdpost;

		LogPoissonInd ();
		virtual ~LogPoissonInd () {}
		void update (const unsigned int&, const arma::mat&,
			arma::ivec&, const arma::mat&, const ParPoissonInd&,
			const PriorPoissonInd&);
};
#endif
