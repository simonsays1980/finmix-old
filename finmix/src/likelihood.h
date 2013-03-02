/* 
 * methods for computation of mixture likelihoods
 * maximum likelihoods and log likelihoods 
 *  author: Lars Simon Zehnder
 * package: finmix
 * created: 02/01/2013
 */

#include <RcppArmadillo.h> 
#include "liklist.h" 		// structure liklist
#include <algorithm> 		// C++ Standard Library algorithms (math functions)
#include <R.h>       		// to interface with R 
#include <Rmath.h>   		// for using internal R C-functions 

inline liklist 
likelihood_poisson(const arma::mat &Y, arma::rowvec lambda) {
	
	/* lambda is a row vector */
	const unsigned int N = Y.n_rows;
	const unsigned int K = lambda.n_elem;
	arma::vec lgammaY(N);
	arma::mat loglik(N, K);
	arma::mat lh(N, K);
		
	for(unsigned int k = 0; k < K; ++k) {
		lambda(k) = std::max(lambda(k), 1e-4);
	}
	
	for(unsigned int i = 0; i < N; ++i) {
		lgammaY(i) = R::lgammafn(Y(i, 0) + 1.0);
	}
	arma::mat lgY = arma::repmat(lgammaY, 1, K);
	arma::mat repY = arma::repmat(Y, 1, K);
	
	for(unsigned int k = 0; k < K; ++k) {
		loglik.row(k) = repY.row(k) % arma::log(lambda);
	}
	loglik.each_row() -= lambda;
	loglik.each_col() -= lgammaY; 

	arma::vec maxl = arma::max(loglik, 1);
	for(unsigned int k = 0; k < K; ++k) {
		lh.col(k) = loglik.col(k) - maxl;
	}

  	liklist l_list(lh, maxl, loglik);

	return l_list;
}
