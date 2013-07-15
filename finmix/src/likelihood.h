/* 
 * methods for computation of mixture likelihoods
 * maximum likelihoods and log likelihoods 
 *  author: Lars Simon Zehnder
 * package: finmix
 * created: 02/01/2013
 */
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <RcppArmadillo.h> 
#include <algorithm> 		// C++ Standard Library algorithms (math functions)
#include <R.h>       		// to interface with R 
#include <Rmath.h>   		// for using internal R C-functions 

struct liklist {
	
	const arma::mat lh;
	const arma::vec maxl;
	const arma::mat llh;
	/* ctor */
	liklist(const arma::mat &lh, const arma::vec &maxl, 
		const arma::mat &llh) : lh(lh), maxl(maxl), llh(llh) {}
};
 
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
	arma::rowvec llambda = arma::log(lambda);
	for(unsigned int i = 0; i < N; ++i) {
		loglik.row(i) = repY.row(i) % llambda;
	}
	loglik.each_row() -= lambda;
	loglik.each_col() -= lgammaY;

	arma::vec maxl = arma::max(loglik, 1);
	for(unsigned int k = 0; k < K; ++k) {
		lh.col(k) = arma::exp(loglik.col(k) - maxl);
	}

  	liklist l_list(lh, maxl, loglik);

	return l_list;
}

inline liklist 
likelihood_poisson (const arma::mat &Y, arma::mat lambda) 
{
	/* lambda is a matrix (exposures in data object) */
	const unsigned int N = Y.n_rows;
	const unsigned int K = lambda.n_cols;
	arma::vec lgammaY(N);
	arma::mat loglik(N, K);
	arma::mat lh(N, K);
	//TODO: Check if using umat with lambda < 1e-04 is faster
	for(unsigned int i = 0; i < N; ++i) {	
		for(unsigned int k = 0; k < K; ++k) {
			lambda(i, k) = std::max(lambda(i, k), 1e-4); 
		}
		lgammaY(i) = R::lgammafn(Y(i, 0) + 1.0);
	}
	arma::mat lgY = arma::repmat(lgammaY, 1, K);
	arma::mat repY = arma::repmat(Y, 1, K);
	arma::mat llambda = arma::log(lambda);
	loglik = repY % llambda;
	loglik -= lambda;
	loglik.each_col() -= lgammaY;
	arma::vec maxl = arma::max(loglik, 1);
	for(unsigned int k = 0; k < K; ++k) {
		lh.col(k) = arma::exp(loglik.col(k) - maxl);
	}
	liklist l_list(lh, maxl, loglik);
	return l_list;
}

inline double  
likelihood_gamma (const arma::rowvec& Y, const double& shape,
	const double& rate) 
{
	const unsigned int N = Y.n_elem;
	double lik = 0.0;
	for(unsigned int i = 0; i < N; ++i) {
		lik += shape * std::log(rate) - R::lgammafn(shape)
			- rate * Y(i) + (shape - 1) * std::log(Y(i)); 
	}
	return lik;
}

inline double 
likelihood_ggamma(const arma::rowvec& lambda, 
	const arma::rowvec& shape, const double& rate,
	const arma::rowvec& loc)
{
	const unsigned int K = lambda.n_elem;
	double lik = 0.0;
	for(unsigned int k = 0; k < K; ++k) {
		lik += shape(k) * std::log(rate) - R::lgammafn(shape(k))
			- rate * (lambda(k) - loc(k)) 
			+ (shape(k) - 1) * std::log(lambda(k) - loc(k));  
	}
	return lik;
}
#endif
