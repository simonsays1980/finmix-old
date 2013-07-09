/*
 * distributions needed for MCMC sampling but not 
 * implemented in Rmath.h
 *  author: Lars Simon Zehnder
 * package: finmix
 * created: 02/01/2013
 */
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
#include <algorithm>		// for use of C++ Standard Library math functions
#include <R.h>
#include <Rmath.h>              // for use of R internal C functions

inline
arma::rowvec rdirichlet (const arma::rowvec& dpar) {
	/**
	 * vector is object of type Col as column-wise
	 * orientation is usual in armadillo 
	 */
	const unsigned int K = dpar.n_elem;
	arma::rowvec par_out(K);
	double sum = 0.0;
	GetRNGstate();

	for(unsigned int k = 0; k < K; ++k) {
		par_out(k) = R::rgamma(dpar(k), 1);
		sum += par_out(k); 
	}

	PutRNGstate();	
	par_out = par_out/sum;

	return par_out;
}

inline
arma::rowvec rgammaprod (const arma::rowvec& par_a, const arma::rowvec& par_b) {
	/** 
	 * vector is object of type Col as column-wise 
	 * orientation is usual in armadillo 
	 */
	const unsigned int K = par_a.n_elem;
	arma::rowvec par_out(K);

	GetRNGstate();
	
	for(unsigned int k = 0; k < K; ++k) {
		par_out(k) = R::rgamma(par_a(k), 1);
		par_out(k) = std::max(par_out(k), 1e-10);
		par_out(k) = par_out(k)/par_b(k); 	
	}

	PutRNGstate();

	return par_out;
}
#endif
