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
#include "PriorCondPoissonFix.h"
#include "ParCondPoissonFix.h"

PriorCondPoissonFix::PriorCondPoissonFix () : PriorPoissonFix() {}


PriorCondPoissonFix::PriorCondPoissonFix (const FinmixPrior& prior) :
	PriorPoissonFix(prior),
	coef(Rcpp::as<arma::mat>((SEXP) prior.par["coef"])),
	cond(Rcpp::as<arma::rowvec>((SEXP) prior.par["a"])) {}

void PriorCondPoissonFix::update (const unsigned int& K, const arma::mat& y,
			arma::ivec& S, const ParPoissonFix& par)  
{
	PriorPoissonFix::update(K, y, S, par);
}

void PriorCondPoissonFix::updateHier(const ParPoissonFix& par)
{
	if (HIER) {
		const unsigned int K = par.lambda.n_elem;
		GetRNGstate();
		double gN = g + arma::sum(aStart);
		double GN = 0.0;
		double b = 0.0;
		for(unsigned int k = 0; k < K; ++k) {
			GN = G * std::pow(2, arma::sum(coef.row(k)) - 1)
				+ arma::sum(par.lambda);
			b = R::rgamma(gN, 1/GN);
			bStart(k) = b;
		}
		PutRNGstate();
	}
}

