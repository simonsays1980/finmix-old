#ifndef ADAPTER_H
#define ADAPTER_H

#include <RcppArmadillo.h>
#include "BASE.h"

template <typename Super>
class ADAPTER : public Super, public BASE {
	public:
		ADAPTER () {}
		ADAPTER (const FinmixData&, const FinmixModel&, const 
			FinmixPrior&, const FinmixMCMC&, const Rcpp::S4&);
		virtual void update ();
		virtual void store (const unsigned int&);
};

template <typename Super>
ADAPTER <Super>::ADAPTER (const FinmixData& data, const FinmixModel& model, const
	FinmixPrior& prior, const FinmixMCMC& mcmc, const Rcpp::S4& classS4) :
		Super(data, model, prior, mcmc, classS4), BASE() {}

template <typename Super>
void ADAPTER <Super>::update () 
{
	Super::update();
}

template <typename Super>
void ADAPTER <Super>::store (const unsigned int& m) 
{
	Super::store(m);
}
#endif
