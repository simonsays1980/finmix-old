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

#ifndef FINMIXDATA_H
#define FINMIXDATA_H

#include <RcppArmadillo.h>
#include <string>

class FinmixData {
	public:
		arma::mat y;
		arma::ivec S;
		arma::vec expos;
		arma::vec T;
	   /**
		* finmix 'data' objects arriving in C++ 
        * are always column-wise ordered 
        * therefore slot 'bycolumn' is left out
        */
		std::string dataType;
		unsigned int N;
		unsigned int r;
	
		/* constructor */ 
		FinmixData(const Rcpp::S4&);
};
#endif
