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


/**
* 'exp' and 'T' can be 'NA' if the model
* has no exposures (Poisson) or repetitions (Binomial)
* in this case these attributes are set to an 1 x 1 vector
*
*/

#include "FinmixData.h"

FinmixData::FinmixData (const Rcpp::S4& classS4) :
	y(Rcpp::as<arma::mat>((SEXP) classS4.slot("y"))),
	S(Rcpp::as<arma::ivec>((SEXP) classS4.slot("S"))),
	expos(Rcpp::as<arma::vec>((SEXP) classS4.slot("exp"))),
	T(Rcpp::as<arma::vec>((SEXP) classS4.slot("T"))),
	dataType(Rcpp::as<std::string>((SEXP) classS4.slot("type"))),
	N(Rcpp::as<unsigned int>((SEXP) classS4.slot("N"))),
	r(Rcpp::as<unsigned int>((SEXP) classS4.slot("r"))) {}		
	
