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
#include "FinmixModel.h"

FinmixModel::FinmixModel(const Rcpp::S4& classS4) {

	/**
	 * As the parameters 'par' can have very differing implementations 
	 * they are kept in a general Rcpp::List object which is then 
	 * decomposed in main code.
	 *
	 */ 
	par  = classS4.slot("par");
	weight  = Rcpp::as<arma::mat>(classS4.slot("weight"));

	/** 
	 * 'T' can be 'NA' if the model has no repetitions (Binomial)
	 * in this case this attribute is set to an 1 x 1 vector
	 *
	 */
	T	 = Rcpp::as<arma::ivec>(classS4.slot("T"));
	indicFix = Rcpp::as<bool>(classS4.slot("indicfix"));
	K 	 = Rcpp::as<unsigned int>(classS4.slot("K"));
	r 	 = Rcpp::as<unsigned int>(classS4.slot("r"));	
}

FinmixModel::~FinmixModel(){}
