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
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <RcppArmadillo.h> 
#include <algorithm> 		// C++ Standard Library algorithms (math functions)
#include <R.h>       		// to interface with R 
#include <Rmath.h>   		// for using internal R C-functions 


/**
 * -----------------------------------------------------------
 * liklist
 * @brief   Structure to hold the return values from the 
 *          likelihood computations. 
 * @author  Lars Simon Zehnder
 * -----------------------------------------------------------
 **/
struct liklist {
	
	const arma::mat lh;
	const arma::vec maxl;
	const arma::mat llh;
	/* ctor */
	liklist(const arma::mat &lh, const arma::vec &maxl, 
		const arma::mat &llh) : lh(lh), maxl(maxl), llh(llh) {}
};

// ===========================================================
// Poisson likelihood
// -----------------------------------------------------------

/**
 * -----------------------------------------------------------
 * likelihood_poisson
 * @brief   Computes likelihood for a Poisson model with 
 *          exposures.
 * @par Y       data values, N x 1
 * @par lambda  parameter vector, 1 x K
 * @return  liklist struct with likelihood values
 * @details the likelihood is computed as well as the log-
 *          likelihood and the maximum of the likelihood 
 *          over components.
 * @see liklist
 * @author Lars Simon Zehnder
 * -----------------------------------------------------------
 **/
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

/**
 * -----------------------------------------------------------
 * likelihood_poisson
 * @brief   Computes likelihood for a Poisson model with 
 *          exposures.
 * @note    this is the default function called 
 * @par Y       data values, N x 1
 * @par lambda  parameter matrix,  N x K
 * @return  liklist struct with likelihood values
 * @details the likelihood is computed as well as the log-
 *          likelihood and the maximum of the likelihood 
 *          over components.
 * @see liklist
 * @author Lars Simon Zehnder
 * -----------------------------------------------------------
 **/
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

// ===========================================================
// Binomial likelihood
// -----------------------------------------------------------

/**
 * -----------------------------------------------------------
 * likelihood_binomial
 * @brief   Computes likelihood for a Binomial model with 
 *          exposures.
 * @note    this is the default function called 
 * @par Y       data values, N x 1
 * @par lambda  parameter matrix,  N x K
 * @par T       repetitions, N x 1
 * @return  liklist struct with likelihood values
 * @details the likelihood is computed as well as the log-
 *          likelihood and the maximum of the likelihood 
 *          over components.
 * @see liklist
 * @author Lars Simon Zehnder
 * -----------------------------------------------------------
 **/

inline
liklist likelihood_binomial (const arma::mat& Y, 
        const arma::rowvec p, 
        const arma::vec& T)
{
   const unsigned int N = Y.n_rows; 
   const unsigned int K = p.n_elem;
   arma::vec lgammaY(N);    
   arma::vec lgammaT(N);
   arma::vec lgammaTY(N);
   arma::mat loglik(N, K);
   arma::mat lh(N, K);   
   for (unsigned int i = 0; i < N; ++i) {
       lgammaY(i) = R::lgammafn(Y(i, 0) + 1.0);      
       lgammaT(i) = R::lgammafn(T(i) + 1.0);
       lgammaTY(i) = R::lgammafn(T(i) - Y(i, 0) + 1.0);
   }   
   arma::mat repY   = arma::repmat(Y, 1, K);
   arma::mat repTY  = arma::repmat(T - Y, 1, K);
   arma::mat lgY    = arma::repmat(lgammaY, 1, K);
   arma::mat lgT    = arma::repmat(lgammaT, 1, K);
   arma::mat lgTY   = arma::repmat(lgammaTY, 1, K);
   for (unsigned int i = 0; i < N; ++i) {
       loglik.row(i) = repY.row(i) % p + repTY.row(i) % arma::log(1.0 - p);
   }
   loglik.each_col() += lgammaT;
   loglik.each_col() -= lgammaTY;
   loglik.each_col() -= lgammaY;
   arma::vec maxl = arma::max(loglik, 1);
   for (unsigned int k = 0; k < K; ++k) {
       lh.col(k) = arma::exp(loglik.col(k) - maxl);
   }
   liklist l_list(lh, maxl, loglik);
   return l_list;
}
#endif
