#ifndef POSTOUTCONDPOISSONIND_H
#define POSTOUTCONDPOISSONIND_H

#include <RcppArmadillo.h>
#include "PostOutCondPoissonFix.h"
#include "PriorCondPoissonInd.h"

class PostOutCondPoissonInd : public PostOutCondPoissonFix {
	public:
		arma::mat* weight;
		
		PostOutCondPoissonInd (const Rcpp::List&);
		virtual ~PostOutCondPoissonInd () {}
		virtual void store (const unsigned int&,
				const PriorCondPoissonInd&); 
};
#endif
