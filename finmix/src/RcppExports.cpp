// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// swap_cc
Rcpp::NumericMatrix swap_cc(Rcpp::NumericMatrix values, Rcpp::IntegerMatrix index);
RcppExport SEXP finmix_swap_cc(SEXP valuesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::NumericMatrix values = Rcpp::as<Rcpp::NumericMatrix >(valuesSEXP);
        Rcpp::IntegerMatrix index = Rcpp::as<Rcpp::IntegerMatrix >(indexSEXP);
        Rcpp::NumericMatrix __result = swap_cc(values, index);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// swapInteger_cc
Rcpp::IntegerMatrix swapInteger_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index);
RcppExport SEXP finmix_swapInteger_cc(SEXP valuesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::IntegerMatrix values = Rcpp::as<Rcpp::IntegerMatrix >(valuesSEXP);
        Rcpp::IntegerMatrix index = Rcpp::as<Rcpp::IntegerMatrix >(indexSEXP);
        Rcpp::IntegerMatrix __result = swapInteger_cc(values, index);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// swapInd_cc
Rcpp::IntegerMatrix swapInd_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index);
RcppExport SEXP finmix_swapInd_cc(SEXP valuesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::IntegerMatrix values = Rcpp::as<Rcpp::IntegerMatrix >(valuesSEXP);
        Rcpp::IntegerMatrix index = Rcpp::as<Rcpp::IntegerMatrix >(indexSEXP);
        Rcpp::IntegerMatrix __result = swapInd_cc(values, index);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// swapST_cc
Rcpp::IntegerVector swapST_cc(Rcpp::IntegerVector values, Rcpp::IntegerMatrix index);
RcppExport SEXP finmix_swapST_cc(SEXP valuesSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::IntegerVector values = Rcpp::as<Rcpp::IntegerVector >(valuesSEXP);
        Rcpp::IntegerMatrix index = Rcpp::as<Rcpp::IntegerMatrix >(indexSEXP);
        Rcpp::IntegerVector __result = swapST_cc(values, index);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
