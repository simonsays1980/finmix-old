/*
 * Rcpp_data.h
 * 
 * definition of Rcppdata struct to contain 
 * finmix 'data' S4-object 
 * Also a converter from Rcpp S4-object to
 * C++-struct is provided.
 *
 *  author: Lars Simon Zehnder
 * package: finmix (1.0.0)
 * created: 02/02/2013
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <string>

template <class M>
class FinmixData {
	public:
		M yM;
		arma::uvec sV;
		arma::vec expV;
		arma::uvec tV;
	       /**
		* finmix 'data' objects arriving in C++ 
         	* are always column-wise ordered 
         	* therefore slot 'bycolumn' is left out
		*
        	*/
		std::string dataType;
		unsigned int N;
		unsigned int R;
	
		/* constructor */ 
		FinmixData(const Rcpp::S4& classS4) {
			/* convert Rcpp objects to Armadillo (C++) objects */
		        dataType = Rcpp::as<std::string>(classS4.slot("type"));

	        	yM 	 = Rcpp::as<M>(classS4.slot("y"));
			sV 	 = Rcpp::as<arma::uvec>(classS4.slot("S"));

        		/** 
		         * 'exp' and 'T' can be 'NA' if the model 
		         * has no exposures (Poisson) or repetitions (Binomial)
        		 * in this case these attributes are set to an 1 x 1 vector
		         *
		         */
		        expV 	 = Rcpp::as<arma::vec>(classS4.slot("exp"));
		        tV 	 = Rcpp::as<arma::uvec>(classS4.slot("T"));

		        N 	 = Rcpp::as<unsigned int>(classS4.slot("N"));
		        R 	 = Rcpp::as<unsigned int>(classS4.slot("r"));
		};	
};

