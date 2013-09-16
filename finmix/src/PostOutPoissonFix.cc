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
#include "PostOutPoissonFix.h"

PostOutPoissonFix::PostOutPoissonFix (const Rcpp::List& list) 
{
	Rcpp::List tmpList((SEXP) list["par"]);
	Rcpp::NumericMatrix tmpA((SEXP) tmpList["a"]);
	Rcpp::NumericMatrix tmpB((SEXP) tmpList["b"]);
	const unsigned int M = tmpA.nrow();
	const unsigned int K = tmpA.ncol();	
	a = new arma::mat(tmpA.begin(), M, K, false, true);
	b = new arma::mat(tmpB.begin(), M, K, false, true);
}

void PostOutPoissonFix::store (const unsigned int& m,
				const PriorPoissonFix& hyperPar) 
{
	(*a).row(m) = hyperPar.aPost;
	(*b).row(m) = hyperPar.bPost;
}
