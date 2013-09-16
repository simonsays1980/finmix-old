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
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
#include <algorithm>		// for use of C++ Standard Library math functions
#include <R.h>
#include <Rmath.h>              // for use of R internal C functions

inline
arma::rowvec rdirichlet (const arma::rowvec& dpar) 
{
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
arma::vec lddirichlet (const arma::mat &values, const arma::vec &par)
{
    const unsigned int M    = values.n_rows;
    const unsigned int K    = values.n_cols;
    arma::vec rvalues       = arma::zeros(M);
    double std_const        = 0.0;
    for (unsigned int k = 0; k < K; ++k) {
        rvalues     += arma::log(values.unsafe_col(k)) * (par.at(k) - 1);
        std_const   += R::lgammafn(par.at(k));
    }
    std_const   -= R::lgammafn(arma::as_scalar(arma::sum(par)));
    rvalues     -= std_const;
    return(rvalues);
}

inline 
arma::rowvec dpoisson(const double &value, const arma::rowvec &par)
{
    const unsigned int K = par.n_elem;
    arma::rowvec rvec(K);
    for (unsigned int k = 0; k < K; ++k) {
        rvec(k) = R::dpois(value, par(k), 0);
    }
    return rvec;
}

inline
arma::rowvec rgammaprod (const arma::rowvec& par_a, 
	const arma::rowvec& par_b) 
{	
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

inline 
arma::mat ldgamma (const arma::mat &values, const arma::vec &shape,
        const arma::vec &rate)
{
    const unsigned int M    = values.n_rows;
    const unsigned int K    = values.n_cols;
    arma::mat rvalues(M, K);
    for (unsigned int k = 0; k < K; ++k) {
        rvalues.unsafe_col(k) = arma::log(values.unsafe_col(k)) * (shape.at(k) - 1);
        rvalues.unsafe_col(k) -= values.unsafe_col(k) * rate.at(k);
        rvalues.unsafe_col(k) += shape.at(k) * std::log(rate.at(k));
        rvalues.unsafe_col(k) -= R::lgammafn(shape.at(k));
    }
    return rvalues;
}

inline 
double rggamma (const double& shape, const double& rate, 
	const double& loc) 
{
	double par_out = 0.0;
	GetRNGstate();
	par_out = R::rgamma(shape, 1);
	PutRNGstate();
	par_out = std::max(par_out, 1e-10);
	par_out = par_out/rate;
	par_out += loc;

	return par_out;
}
#endif
