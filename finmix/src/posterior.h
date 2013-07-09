/**
 * computes the moments of a posterior distribution under a 
 * conditionally conjugate prior
 *  author: Lars Simon Zehnder
 * package: 'finmix'
 * created: 31 Jan. 2013
 */
#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <RcppArmadillo.h>

/**
 * posterior for multinomial distribution
 *
 * as this is used for the weights of any mixture 
 * no specific struct is input argument 
 */

inline arma::rowvec 
posterior_multinomial(const unsigned int &K, const arma::ivec &S, 
			const arma::rowvec &weight) 
{
	arma::imat repS = arma::repmat(S, 1, K);
	arma::imat compM = arma::ones<arma::imat>(S.n_rows, K);
	arma::rowvec par_post(K);
 
	/* create sequence */
	for(unsigned int k = 0; k < K; ++k) {
		compM.col(k) = compM.col(k) * (k + 1);		
	}
	arma::umat ind = (repS == compM);
	arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
	par_post = sum(indDouble); 
	par_post = par_post + weight;
	
	return par_post;
}
#endif
