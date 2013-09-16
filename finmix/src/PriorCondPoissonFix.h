/******************************************************************************
 *
 * Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is part of the R package finmix.
 *
 * finmix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundatio, either version 3 of the License, or
 * any later version.
 *
 * finmix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with finmix. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef PRIORCONDPOISSONFIX_H
#define PRIORCONDPOISSONFIX_H

#include <RcppArmadillo.h>
#include "PriorPoissonFix.h"
#include "FinmixPrior.h"

class PriorCondPoissonFix : virtual public PriorPoissonFix {
	public:
		arma::mat coef;
		arma::rowvec cond;

		PriorCondPoissonFix ();	
		PriorCondPoissonFix	(const FinmixPrior&);
		virtual	~PriorCondPoissonFix () {} 
		virtual void update (const unsigned int&, 
			const arma::mat&, arma::ivec&,
			const ParPoissonFix&);
		virtual void updateHier(const ParPoissonFix&);
};
#endif
