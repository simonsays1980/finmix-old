/******************************************************************************
 *
 * Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is part of the R package 'finmix'.
 *
 * 'finmix' is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundatio, either version 3 of the License, or
 * any later version.
 *
 * 'finmix' is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/


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
