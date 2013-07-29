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
#ifndef MCMCCONDPOISSON_CC
#define MCMCCONDPOISSON_CC

#include <RcppArmadillo.h>		// C++ linear algebra library
#include "FinmixData.h"
#include "FinmixModel.h"
#include "FinmixPrior.h"
#include "FinmixMCMC.h"
#include "BASE.h"
#include "ADAPTER.h"
#include "FIX.h"
#include "IND.h"
#include "HIER.h"
#include "POST.h"
#include "LogCondPoissonFix.h"
#include "LogCondPoissonInd.h"
#include "ParCondPoissonInd.h"
#include "ParOutPoisson.h"
#include "HierOutPoisson.h"
#include "PostOutCondPoissonFix.h"
#include "PostOutCondPoissonInd.h"


RcppExport SEXP mcmc_condpoisson_cc(SEXP data_S4, SEXP model_S4, 
	SEXP prior_S4, SEXP mcmc_S4, SEXP mcmcoutput_S4) 
{
	
	/* Convert S4-classes to C++-structs */
	Rcpp::S4 dataS4O(data_S4);
	Rcpp::S4 modelS4O(model_S4);
	Rcpp::S4 priorS4O(prior_S4);
	Rcpp::S4 mcmcS4O(mcmc_S4);
	Rcpp::S4 mcmcOutputS4O(mcmcoutput_S4);
	FinmixData finData = FinmixData(dataS4O);
	FinmixModel finModel = FinmixModel(modelS4O);
	FinmixPrior finPrior = FinmixPrior(priorS4O);
	FinmixMCMC finMCMC = FinmixMCMC(mcmcS4O);

	const bool INDICFIX 		= finModel.indicFix;
	const bool HIER_IND 		= finPrior.hier;
	const bool POST_IND 		= finMCMC.storePost;
	const unsigned int BURNIN 	= finMCMC.burnIn;
	const unsigned int M 		= finMCMC.M;
	const unsigned int K 		= finModel.K;
	
	BASE* ptr;
	typedef FIX<PriorCondPoissonFix, ParCondPoissonFix, 
		LogCondPoissonFix, ParOutPoisson> CONDPOISSONFIX;
	typedef IND<FIX<PriorCondPoissonInd, ParCondPoissonInd, 
		LogCondPoissonInd, ParOutPoisson> > CONDPOISSONIND;
	if (INDICFIX || K == 1) 
	{ 
		if (HIER_IND) {
			if (POST_IND) {
					ptr = new ADAPTER<POST<HIER<CONDPOISSONFIX,
					HierOutPoisson>, PostOutCondPoissonFix> >
					(finData, finModel, finPrior, finMCMC,
					mcmcOutputS4O);
	
			} else {
				ptr = new ADAPTER<HIER<CONDPOISSONFIX, HierOutPoisson> >
					(finData, finModel, finPrior, finMCMC, 
					mcmcOutputS4O);
			}
		} else {
			if (POST_IND) {
				ptr = new ADAPTER<POST<CONDPOISSONFIX,
						PostOutCondPoissonFix> >
					(finData, finModel, finPrior, finMCMC, 
					mcmcOutputS4O);
			} else {
				ptr = new ADAPTER<CONDPOISSONFIX> (finData, finModel,
					finPrior, finMCMC, mcmcOutputS4O);
			}
		}
	} else {	
		if (HIER_IND) {
			if (POST_IND) {
				ptr = new ADAPTER<POST<HIER<CONDPOISSONIND,
					HierOutPoisson>, PostOutCondPoissonInd> >
					(finData, finModel, finPrior, finMCMC,
					mcmcOutputS4O);
			} else {
				ptr = new ADAPTER<HIER<CONDPOISSONIND, HierOutPoisson> >
					(finData, finModel, finPrior, finMCMC,
					mcmcOutputS4O);
			} 
		} else {
			if (POST_IND) {
				ptr = new ADAPTER<POST<CONDPOISSONIND, PostOutCondPoissonInd> >
					(finData, finModel, finPrior, finMCMC,
					mcmcOutputS4O);
			} else {
				ptr = new ADAPTER<CONDPOISSONIND> (finData, finModel,
					finPrior, finMCMC, mcmcOutputS4O);
			}
		}
	}
	for(unsigned int i = 0; i < BURNIN + M; ++i) {
		ptr->update();
		ptr->store(i);
	}		
		
	return Rcpp::wrap(mcmcOutputS4O);	
}
#endif
