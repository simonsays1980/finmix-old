#ifndef LOGCONDPOISSONIND_H
#define LOGCONDPOISSONIND_H

#include <RcppArmadillo.h>
#include "LogCondPoissonFix.h"
#include "ParCondPoissonInd.h"
#include "PriorCondPoissonInd.h"

class LogCondPoissonInd : public LogCondPoissonFix {
	public:
		double cdpost;
		double entropy;
		double maxcdpost;

		LogCondPoissonInd ();
		virtual ~LogCondPoissonInd () {}
		void update (const unsigned int&, 
			const arma::mat&, arma::ivec&, 
			const arma::mat&, 
			const ParCondPoissonInd&,
			const PriorCondPoissonInd&);
};
#endif
