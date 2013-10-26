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
#ifndef __FINMIX_DISTRIBUTIONS_H__
#define __FINMIX_DISTRIBUTIONS_H__

#include <RcppArmadillo.h>
#include <algorithm>		// for use of C++ Standard Library math functions
#include <R.h>
#include <Rmath.h>              // for use of R internal C functions

// ===============================================================
// Dirichlet distribution
// ---------------------------------------------------------------

/**
 * ---------------------------------------------------------------
 * @brief   Samples a vector from a Dirichlet distribution.
 * @par dpar    Dirichlet parameters, 1 x K
 * @detail  Sampling of the Dirichlet is implemented through
 *          the function 'Rf_rgamma()' from Rmath.h. 
 * @see R::rgamma
 * @author Lars Simon Zehnder
 * ---------------------------------------------------------------
 **/
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

/**
 * ------------------------------------------------------------------
 * lddirichlet
 * @brief   Computes the log density of the Dirichlet distribution 
 *          for a vector of values.
 * @param   values  values for which the log density should be 
 *          calculated
 * @param   par     parameters of the Dirichlet distribution
 * @detail  The function does use the fast access function 'at()' 
 *          for Armadillo objects and therefore no boundaries for
 *          indices get checked. Inside the Rcpp sugar wrapper 
 *          'lgammafn()' from the 'R' namespace is used. 
 * @see R::lgammfn, arma::mat<>.at()
 * @author Lars Simon Zehnder
 * -----------------------------------------------------------------
 */
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

// =================================================================
// Poisson distribution
// -----------------------------------------------------------------

/**
 * -----------------------------------------------------------------
 * dpoisson
 * @brief   Computes density of the Poisson distribution for an
 *          Armadillo parameter vector.
 * @param   value   the density is calculated for
 * @param   par     parameter vector; 1 x K
 * @return  vector with density values for the corresponding 
 *          parameters in 'par'
 * @detail  Uses inside the 'dpois()' function from Rcpp's 'R'
 *          namespace.
 * @see R::dpois
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------------
 **/
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

// ==============================================================
// Gamma distribution
// --------------------------------------------------------------

/** 
 * --------------------------------------------------------------
 * @brief   Draws vector random sample for Gamma distribution.
 * @par par_a   arma::vec with shape parameters
 * @par par_b   arma::vec with rate parameters
 * @return  random sample from Gamma distribution with para-
 *          meters defined in par_a and par_b.
 * @detail  uses the C rgamma() function of R
 * @see     R::rgamma
 * @author Lars Simon Zehnder
 * --------------------------------------------------------------
 **/
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

/**
 * ---------------------------------------------------------------
 * ldgamma
 * @brief   Computes the log density of the Gamma distri-
 *          bution for a vector of values.
 * @par values  values for which the log density should be cal-
 *              culated; M x K
 * @par shape   Gamma shape parameters; K x 1
 * @par rate    Gamma rate parameters; K x 1
 * @return  Armadillo matrix with the log densities for each value
 *          in a row and for each pair of parameters in a column;
 *          M x K
 * @detail  For each shape and rate parameter pair the log gamma 
 *          density is computed. Inside the function the unsafe
 *          access functions of Armadillo 'at()' and 'unsafe_col()'
 *          are used, so now boundary check is performed. In each
 *          step the 'lngamma()' function from Rcpp's 'R' namespace
 *          is used.
 * @see R::lgammafn, arma::vec<>::at(), arma::mat<>::unsafe_col
 * @author  Lars Simon Zehnder
 * ----------------------------------------------------------------
 **/
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

// =======================================================
// Beta distribution
// -------------------------------------------------------

/**
 * -------------------------------------------------------
 * rbetaprod
 * @brief   Draws a random vector sample from the Beta 
 *          distribution.
 * @par par_a   shape parameters
 * @par par_b   shape parameters
 * @return  random sample from Beta distribution with para-
 *          meters defined in par_a and par_b
 * @detail  uses the C rbeta() function of R
 * @see     Rf_rbeta
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
inline
arma::rowvec rbetaprod (const arma::rowvec& par_a, 
        const arma::rowvec& par_b)
{
    const unsigned int K = par_a.n_elem;
    arma::rowvec par_out(K);
    Rcpp::RNGScope scope;
    for (unsigned int k = 0; k < K; ++k) {
        par_out(k) = R::rbeta(par_a(k), par_b(k));
        par_out(k) = std::max(par_out(k), 1e-10);
    }
    return par_out;
}

/**
 * ---------------------------------------------------------------
 * ldbeta
 * @brief   Computes the log density of the Beta distri-
 *          bution for a vector of values.
 * @par values  values for which the log density should be cal-
 *              culated; M x K
 * @par shape   Beta first shape parameters; K x 1
 * @par rate    Beta second parameters; K x 1
 * @return  Armadillo matrix with the log densities for each value
 *          in a row and for each pair of parameters in a column;
 *          M x K
 * @detail  For each shape1 and shape2 parameter pair the log Beta 
 *          density is computed. Inside the function the unsafe
 *          access functions of Armadillo 'at()' and 'unsafe_col()'
 *          are used, so now boundary check is performed. In each
 *          step the 'lbeta()' function from Rcpp's 'R' namespace
 *          is used.
 * @see R::lbeta, arma::vec<>::at(), arma::mat<>::unsafe_col
 * @author  Lars Simon Zehnder
 * ----------------------------------------------------------------
 **/

inline 
arma::mat ldbeta (const arma::mat &values, const arma::vec &shape1,
        const arma::vec &shape2)
{
    const unsigned int M    = values.n_rows;
    const unsigned int K    = values.n_cols;
    arma::mat rvalues(M, K);
    for (unsigned int k = 0; k < K; ++k) {
        rvalues.unsafe_col(k) = arma::log(values.unsafe_col(k)) * (shape1.at(k) - 1);
        rvalues.unsafe_col(k) += arma::log(values.unsafe_col(k)) * (shape2.at(k) - 1);
        rvalues.unsafe_col(k) -= R::lbeta(shape1.at(k), shape2.at(k));
    }
    return rvalues;
}

// =================================================================
// Binomial distribution
// -----------------------------------------------------------------

/**
 * -----------------------------------------------------------------
 * dpoisson
 * @brief   Computes density of the Binomial distribution for an
 *          Armadillo parameter vector.
 * @param   value   the density is calculated for
 * @param   T       repetitions for the Binomial distribution
 * @param   par     parameter vector; 1 x K
 * @return  vector with density values for the corresponding 
 *          parameters in 'par'
 * @detail  Uses inside the 'dpois()' function from Rcpp's 'R'
 *          namespace.
 * @see R::dpois
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------------
 **/
inline 
arma::rowvec dbinomial(const double& value, const double& T, 
        const arma::rowvec& par)
{
    const unsigned int K = par.n_elem;
    arma::rowvec rvec(K);
    for (unsigned int k = 0; k < K; ++k) {
        rvec(k) = R::dbinom(value, T, par(k), 0);
    }
    return rvec;
}

#endif // __FINMIX_DISTRIBUTIONS_H__
