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
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "distributions.h"
#include "algorithms.h"
#include "optimize.h"
#include "hungarian.h"
#include <nlopt.hpp>
// [[Rcpp::export]]

Rcpp::NumericMatrix swap_cc(Rcpp::NumericMatrix values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if(values.nrow() != index.nrow() || values.ncol() != index.ncol()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Do not reuse memory from R as otherwise existing objects 
     * get manipulated */
    const unsigned int K = values.ncol();
    const unsigned int M = values.nrow();
    arma::mat values_arma(values.begin(), M, K, false, true);
    arma::imat index_arma(index.begin(), M, K, false, true);
    arma::mat values_copy(M, K);
    arma::umat index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
    arma::uvec row_index(1);
    arma::urowvec swap_index(K);
    for(unsigned int i = 0; i < M; ++i) {
        row_index.at(0) = i;
        swap_index = index_umat.row(i);
        values_copy.row(i) = 
            values_arma.submat(row_index, swap_index);
    }
    return Rcpp::wrap(values_copy);
}


// [[Rcpp::export]]

Rcpp::IntegerMatrix swapInteger_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if(values.nrow() != index.nrow() || values.ncol() != index.ncol()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Do not reuse memory from R as otherwise existing objects 
     * get manipulated */
    const unsigned int K = values.ncol();
    const unsigned int M = values.nrow();
    arma::imat values_arma(values.begin(), M, K, false, true);
    arma::imat index_arma(index.begin(), M, K, false, true);
    arma::imat values_copy(M, K);
    arma::umat index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
    arma::uvec row_index(1);
    arma::urowvec swap_index(K);
    for(unsigned int i = 0; i < M; ++i) {
        row_index.at(0) = i;
        swap_index = index_umat.row(i);
        values_copy.row(i) = 
            values_arma.submat(row_index, swap_index);
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]

Rcpp::IntegerMatrix swapInd_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if (values.ncol() != index.nrow()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Reuse memory from R */
    const unsigned int N = values.nrow();
    const unsigned int STORES = values.ncol();
    const unsigned int M = index.nrow();
    const unsigned int K = index.ncol();
    arma::imat values_arma(values.begin(), N, STORES, false, true);
    arma::imat index_arma(index.begin(), M, K, false, true);
    arma::imat values_copy(N, STORES);
    for (unsigned int s = 0; s < STORES; ++s) {
        for (unsigned int i = 0; i < N; ++i) {
            values_copy(i, s) = (int) index_arma(s, (unsigned int)
                    values_arma(i, s) - 1);
        }
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]

Rcpp::IntegerVector swapST_cc(Rcpp::IntegerVector values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if (values.size() != index.nrow()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Reuse memory from R */
    const unsigned int M = values.size();
    const unsigned int K = index.ncol();
    arma::ivec values_arma(values.begin(), M, false, true);
    arma::imat index_arma(index.begin(), M, K, false, true);
    arma::ivec values_copy(M);
    for(unsigned int i = 0; i < M; ++i) {
        values_copy(i) = index_arma(i, (unsigned int)
                values_arma(i) - 1);
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]

Rcpp::NumericMatrix ldgamma_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector shape,
        Rcpp::NumericVector rate)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_shape(shape.begin(), K, false, true);
    arma::vec arma_rate(rate.begin(), K, false, true);
    arma::mat arma_return(M, K);
    arma_return = ldgamma(arma_values, arma_shape, arma_rate);
    return Rcpp::wrap(arma_return);
}

// [[Rcpp::export]]

arma::mat dgamma_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector shape,
        Rcpp::NumericVector rate)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_shape(shape.begin(), K, false, true);
    arma::vec arma_rate(rate.begin(), K, false, true);
    arma::mat arma_return(M, K);
    arma_return = exp(ldgamma(arma_values, arma_shape, arma_rate));
    return arma_return;
}

// [[Rcpp::export]]

Rcpp::NumericVector lddirichlet_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector par)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_par(par.begin(), K, false, true);
    arma::vec arma_return(M);
    arma_return = lddirichlet(arma_values, arma_par);
    return Rcpp::wrap(arma_return);
}

// [[Rcpp::export]]

arma::vec ddirichlet_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector par)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_par(par.begin(), K, false, true);
    arma::vec arma_return(M);
    arma_return = arma::exp(lddirichlet(arma_values, arma_par));
    return arma_return;
}

// [[Rcpp::export]]

arma::imat maxlabel_poisson_cc(const arma::mat values1, const arma::mat values2, 
        const arma::vec shape, const arma::vec rate, const arma::vec dirich,
        const arma::umat perm)
{
    unsigned int M  = values1.n_rows;
    unsigned int K  = values1.n_cols;
    unsigned int P  = perm.n_rows;
    arma::umat index(M, K);
    arma::mat func_val(M, P);
    arma::uvec row_index(M);
    arma::uvec col_index(P);
    arma::vec tmp(M);
    arma::umat arma_perm = perm - 1;
    for (unsigned int m = 0; m < M; ++m) {
        row_index(m) = m;
    }
    for (unsigned int p = 0; p < P; ++p) {                   
        tmp = arma::prod(arma::exp(ldgamma(values1(row_index, 
                            arma_perm.row(p)), shape, rate)), 1);
        tmp %= arma::exp(lddirichlet(values2(row_index, arma_perm.row(p)), 
                    dirich));
        func_val.col(p) = arma::log(tmp);
    }
    for (unsigned int m = 0; m < M; ++m) {    
        arma::vec tmp   = arma::conv_to<arma::vec>::from(func_val.row(m));
        col_index       = arma::sort_index(tmp, 1);
        index.row(m)    = arma_perm.row(col_index(0));
    }
    index = index + 1;
    return arma::conv_to<arma::imat>::from(index);
}

// [[Rcpp::export]] 

arma::imat hungarian_cc(const arma::mat cost) 
{
    arma::umat indM = hungarian(cost);
    return arma::conv_to<arma::imat>::from(indM);
}

// [[Rcpp::export]]

arma::imat stephens1997a_poisson_cc(const Rcpp::NumericMatrix values1, 
        const Rcpp::NumericMatrix values2, 
        arma::vec pars, const arma::umat perm)
{
    const unsigned int M = values1.rows();
    const unsigned int K = values2.cols();
    const unsigned int P = perm.n_rows;
    const unsigned int n = pars.n_elem;
    double value = 1.0;
    double value_next = 0.0;
    arma::mat lambda(values1.begin(), M, K, true, true);
    arma::mat weight(values2.begin(), M, K, true, true);
    const arma::umat arma_perm = perm - 1;
    arma::uvec row_index(M);
    arma::uvec col_index(K);
    arma::vec tmp(M);
    arma::vec tmp2(K);
    arma::umat index = arma::ones<arma::umat>(M, K);
    arma::umat ind(M, K);
    arma::vec dirich(K);
    arma::vec shape(K);
    arma::vec rate(K);
    arma::mat func_val(M, K);
    for (unsigned int k = 0; k < K; ++k) {
        index.unsafe_col(k) *= k;
    }
    for (unsigned int m = 0; m < M; ++m) {
        row_index.at(m) = m;
    }
    /* Set up the optimizer */
    nlopt::opt optim(nlopt::LN_NELDERMEAD, n);
    std::vector<arma::mat*> f_data(2);
    f_data[0] = &lambda;
    f_data[1] = &weight;
    optim.set_max_objective(obj_stephens1997a_poisson, &f_data);
    optim.set_lower_bounds(1e-10);
    optim.set_upper_bounds(1e+7);
    std::vector<double> opt_par = arma::conv_to<std::vector<double> >::from(pars);

    while (value != value_next) {
        value = value_next;
        nlopt::result result = optim.optimize(opt_par, value_next);
        for (unsigned int k = 0; k < K; ++k) {
            dirich.at(k)   = opt_par[k];
            shape.at(k)    = opt_par[k + K];
            rate.at(k)     = opt_par[k + 2 * K];
        }
        /* Loop over permutations */
        for (unsigned int p = 0; p < P; ++p) {
            tmp = arma::prod(arma::exp(ldgamma(lambda(row_index, arma_perm.row(p)), shape, rate)), 1);
            tmp %= arma::exp(lddirichlet(weight(row_index, arma_perm.row(p)), dirich));
            func_val.unsafe_col(p) = arma::log(tmp);
        }
        for (unsigned int m = 0; m < M; ++m) {
            tmp2            = arma::conv_to<arma::vec>::from(func_val.row(m));
            col_index       = arma::sort_index(tmp2, 1);
            ind.row(m)    = arma_perm.row(col_index(0));
        }
        swapmat_by_index(lambda, ind);
        swapmat_by_index(weight, ind);
        swapumat_by_index(index, ind);
    }
    index += 1;
    return arma::conv_to<arma::imat>::from(index);
}

// [[Rcpp::export]]

arma::imat stephens1997b_poisson_cc(const Rcpp::NumericVector values, 
        const Rcpp::NumericMatrix comp_par, 
        const Rcpp::NumericMatrix weight_par)
{
    unsigned int N      = values.size();
    unsigned int M      = comp_par.rows();
    unsigned int K      = comp_par.cols();
    double value = 1.0;
    double value_next = 0.0;
    arma::vec arma_values(values.begin(), N, false, true); 
    arma::mat lambda(comp_par.begin(), M, K, true, true);
    arma::mat weight(weight_par.begin(), M, K, true, true);
    arma::umat index(M, K);
    arma::umat index_out(M, K);
    arma::umat indM(K, K); 
    arma::mat pmat_hat(N, K);
    arma::mat cost(K, K);
    arma::uvec seq_vec(K);
    std::vector<arma::mat*> mat_vector(M);
    pmat_hat    = arma::zeros(N, K);
    index_out   = arma::ones<arma::umat>(M, K);
    for (unsigned int k = 0; k < K; ++k) {
        seq_vec.at(k) = k * K;
        index_out.unsafe_col(k) *= (k + 1);
    }
    for (unsigned int m = 0; m < M; ++m) {
        arma::mat* pmat_ptr = new arma::mat(N, K);
        /* Save a pointer to the STL vector */
        mat_vector[m] = pmat_ptr;
    }
    while (value != value_next) {
        value       = value_next;
        value_next  = 0.0;
        /* For all sampled MCMC parameters a matrix 
         * pmat (N x K) is computed with p_ij 
         * indicating the probability for a value i
         * being from component j.
         * */
        for (unsigned int m = 0; m < M; ++m) {
            for (unsigned int n = 0; n < N; ++n) {
                mat_vector[m]->row(n) = weight.row(m) 
                    % dpoisson(arma_values.at(n), lambda.row(m)); 
                mat_vector[m]->row(n) /= arma::sum(mat_vector[m]->row(n));
            }
        }
        for (unsigned int m = 0; m < M; ++m) {
            pmat_hat += *(mat_vector[m]);
        }
        /* This computes the reference estimator P_hat*/
        pmat_hat /= M;
        /* Now for each sampled MCMC parameter it is
         * searched for the optimal label by computing 
         * the Kullback-Leibler distance of each 'pmat'
         * column 'l' from column 'k' of the reference 
         * estimator P_hat.
         * The cost matrix cost_mat contains then the 
         * distance of column 'l' from column 'k'. 
         * An optimal assignment method computes the minimal
         * 'cost' regarding the labeling.
         * If 'k' is therein assigned to 'l', than the
         * label 'k' is switched to 'l'.
         * */
        for (unsigned int m = 0; m < M; ++m) {
            for (unsigned int k = 0; k < K; ++k) {
                for (unsigned int l = 0; l < K; ++l) {                
                    arma::vec mycol = mat_vector[m]->unsafe_col(l);
                    cost(k, l) = kulback_leibler(mat_vector[m]->unsafe_col(l), 
                            pmat_hat.unsafe_col(k));
                }          
            }
            value_next += arma::trace(cost);
            /* Assignment */
            indM = hungarian(cost);
            arma::uvec f = arma::find(indM.t() == 1);
            index.row(m) = arma::trans(arma::find(indM.t() == 1) - seq_vec);  
        }        
        /* Permute parameters */
        swapmat_by_index(lambda, index);
        swapmat_by_index(weight, index);        
        swapumat_by_index(index_out, index);
        pmat_hat = arma::zeros(N, K);
    }        
    for (unsigned int m = 0; m < M; ++m) {
        delete mat_vector[m];
    }
    return arma::conv_to<arma::imat>::from(index_out);
}
