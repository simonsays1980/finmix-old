/******************************************************************************
 *
 * TODO: Project Title
 *
 * Copyright (C) 2012-2013 Lars Simon Zehnder. All Rights Reserved.
 * Web: -
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 ******************************************************************************/

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export]]

Rcpp::NumericMatrix swap_cc(Rcpp::NumericMatrix values, Rcpp::IntegerMatrix index) {
    /* If both arguments have not the same dimension throw an exception*/
    if(values.nrow() != index.nrow() || values.ncol() != index.ncol()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Reuse memory from R */
    const unsigned int K = values.ncol();
    const unsigned int M = values.nrow();
    arma::mat values_arma(values.begin(), M, K, false, true);
    arma::imat index_imat(index.begin(), M, K, false, true);
    arma::umat index_arma = arma::conv_to<arma::umat>::from(index_imat);
    arma::urowvec comp_index = arma::ones<arma::umat>(K);
    for(unsigned int k = 0; k < K; ++k) {
        comp_index(k) = comp_index(k) * (k + 1); 
    }
    
    /* Swap each row */
    arma::urowvec row_ind(1);
    for(unsigned int i = 0; i < M; ++i) {
        row_ind(1) = i;
        values_arma.submat(row_ind, comp_index) = 
            values_arma.submat(row_ind, index_arma.row(i));
    }
    
    return values;
}

// [[Rcpp::export]]

Rcpp::NumericMatrix swap_cc2(Rcpp::NumericMatrix values, Rcpp::IntegerMatrix index) {
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
    arma::imat compM = arma::ones<arma::imat>(M, K);
    arma::umat ind(M, K);
    arma::mat indDouble(M, K);
    for(unsigned int k = 0; k < K; ++k) {
        ind = (index_arma == (k + 1) * compM);
        indDouble = arma::conv_to<arma::mat>::from(ind);
        values_copy.col(k) = arma::sum(values_arma % indDouble, 1);
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
    arma::urowvec swapIndex(K);
    for (unsigned int s = 0; s > STORES; ++s) {
        swapIndex = arma::sort_index(index_arma.row(s));
        for(unsigned int i = 0; i < N; ++i) {
            values_copy(s, i) = (int) swapIndex((unsigned int)
                    values_arma(s, i) - 1) + 1;
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
                values_arma(i) - 1) + 1;
    }
    return Rcpp::wrap(values);
}
