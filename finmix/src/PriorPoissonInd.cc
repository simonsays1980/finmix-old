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
#include "PriorPoissonInd.h"
#include "ParPoissonInd.h"
#include "posterior.h"

PriorPoissonInd::PriorPoissonInd (const FinmixPrior& prior) :
	PriorPoissonFix(prior),
	weightStart(prior.weight),
	weightPost(prior.weight) {}
void PriorPoissonInd::update (const unsigned int& K, const arma::mat& y,
	arma::ivec& S, const ParPoissonInd& par) 

{
	PriorPoissonFix::update(K, y, S, par);
	weightPost = posterior_multinomial(K, S, weightStart);	
}
