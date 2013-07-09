/**
 * This class holds all the output data
 * 
 * author: Lars Simon Zehnder
 * package: finmix
 * created: 4 March 2013
 *
 */

#include <RcppArmadillo.h>

class FinmixMCMCOutput {
	
	public: 
		unsigned int M;
		Rcpp::NumericMatrix weightsM;
		Rcpp::List par;
		bool ranPerm;
		Rcpp::List hyper;
		Rcpp::List log;
		Rcpp::NumericVector entropyV;
		Rcpp::NumericVector sTM;
		Rcpp::IntegerMatrix sM;
		Rcpp::NumericMatrix nKs;
		Rcpp::List posts;
		Rcpp::IntegerVector clustV;
		
		FinmixMCMCOutput(Rcpp::S4& classS4);
		~FinmixMCMCOutput();
};
