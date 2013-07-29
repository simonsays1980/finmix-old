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

#include "ParCondPoissonFix.h"

ParCondPoissonFix::ParCondPoissonFix (const bool& STARTPAR, 
	const FinmixModel& model) : 
		ParPoissonFix::ParPoissonFix(true, model) {}

void ParCondPoissonFix::update (PriorCondPoissonFix& hyperPar) 
{
	const unsigned int K = lambda.n_elem;
	double tmp = 0.0;
	for(unsigned int k = 0; k < K; ++k) {
		tmp = arma::as_scalar(hyperPar.coef.row(k) 
			* lambda.t());
		tmp -= lambda(k);
		hyperPar.cond(k) = tmp;	
		lambda(k) = rggamma(hyperPar.aPost(k), 
			hyperPar.bPost(k), tmp);
	}
}
