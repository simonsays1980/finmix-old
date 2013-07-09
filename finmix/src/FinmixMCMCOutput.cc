/**
 *
 *
 *
 */

#include <RcppArmadillo.h> 
#include "FinmixMCMCOutput.h"

FinmixMCMCOutput::FinmixMCMCOutput (Rcpp::S4& classS4) : 
			M(Rcpp::as<unsigned int>(classS4.slot("M"))),
			weightsM((SEXP) classS4.slot("weight")),
			par((SEXP) classS4.slot("par")),
			ranPerm(Rcpp::as<bool>(classS4.slot("ranperm"))),
			hyper((SEXP) classS4.slot("hyper")),
			log((SEXP) classS4.slot("log")),
			entropyV((SEXP) classS4.slot("entropy")),
			sTM((SEXP) classS4.slot("ST")),
			sM((SEXP) classS4.slot("S")),
			nKs((SEXP) classS4.slot("NK")),
			posts((SEXP) classS4.slot("post")),
			clustV((SEXP) classS4.slot("clust")){}
			
FinmixMCMCOutput::~FinmixMCMCOutput(){}
