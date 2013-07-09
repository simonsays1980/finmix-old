#ifndef POSTOUTPOISSONIND_H
#define POSTOUTPOISSONIND_H

#include <RcppArmadillo.h>
#include "PostOutPoissonFix.h"
#include "PriorPoissonInd.h"

class PostOutPoissonInd : public PostOutPoissonFix {
	public:
		arma::mat* weight;
		
		PostOutPoissonInd (const Rcpp::List&);
		virtual ~PostOutPoissonInd () {}
		virtual void store (const unsigned int&,
				const PriorPoissonInd&); 
};
#endif
