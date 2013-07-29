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
 * along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
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
