/*
 * structure containing 
 * lh  : likelihood
 * maxl: maximum log likelihood
 * llh : log likelihood
 *
 * author : Lars Simon Zehnder
 * package: finmix
 * created: 02/01/2013
 */

#include <RcppArmadillo.h> 

struct liklist {
	
	const arma::mat lh;
	const arma::vec maxl;
	const arma::mat llh;
	/* ctor */
	liklist(const arma::mat &lh, const arma::vec &maxl, 
		const arma::mat &llh) : lh(lh), maxl(maxl), llh(llh) {}
}; 
