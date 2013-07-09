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
#ifndef FINMIXDATA_H
#define FINMIXDATA_H

#include <RcppArmadillo.h>
#include <string>

class FinmixData {
	public:
		arma::mat y;
		arma::ivec S;
		arma::vec expos;
		arma::ivec T;
	       /**
		* finmix 'data' objects arriving in C++ 
         	* are always column-wise ordered 
         	* therefore slot 'bycolumn' is left out
		*
        	*/
		std::string dataType;
		unsigned int N;
		unsigned int r;
	
		/* constructor */ 
		FinmixData(const Rcpp::S4&);
};
#endif
