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
#ifndef FINMIXMCMC_H
#define FINMIXMCMC_H

#include <RcppArmadillo.h>

class FinmixMCMC {
	
	public:
		unsigned int burnIn;
		unsigned int M;
		bool startPar;
		unsigned int storeS;
		bool storePost;
		bool ranPerm;
	
		/* ctor */ 
		FinmixMCMC(const Rcpp::S4& classS4);	
		/* dtor */
		~FinmixMCMC();
};
#endif
