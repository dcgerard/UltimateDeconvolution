// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dnorm_rank1
double dnorm_rank1(NumericVector x, NumericVector v, NumericVector s_diag, NumericVector mu, bool return_log);
RcppExport SEXP UltimateDeconvolution_dnorm_rank1(SEXP xSEXP, SEXP vSEXP, SEXP s_diagSEXP, SEXP muSEXP, SEXP return_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_diag(s_diagSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type return_log(return_logSEXP);
    rcpp_result_gen = Rcpp::wrap(dnorm_rank1(x, v, s_diag, mu, return_log));
    return rcpp_result_gen;
END_RCPP
}
// get_llike_mat_cpp
NumericMatrix get_llike_mat_cpp(NumericMatrix x_mat, NumericMatrix s_mat, NumericMatrix v_mat, NumericVector pi_vec);
RcppExport SEXP UltimateDeconvolution_get_llike_mat_cpp(SEXP x_matSEXP, SEXP s_matSEXP, SEXP v_matSEXP, SEXP pi_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s_mat(s_matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type v_mat(v_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(get_llike_mat_cpp(x_mat, s_mat, v_mat, pi_vec));
    return rcpp_result_gen;
END_RCPP
}
// dmixlike_cpp
double dmixlike_cpp(NumericMatrix x_mat, NumericMatrix s_mat, NumericMatrix v_mat, NumericVector pi_vec, bool return_log);
RcppExport SEXP UltimateDeconvolution_dmixlike_cpp(SEXP x_matSEXP, SEXP s_matSEXP, SEXP v_matSEXP, SEXP pi_vecSEXP, SEXP return_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s_mat(s_matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type v_mat(v_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_vec(pi_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type return_log(return_logSEXP);
    rcpp_result_gen = Rcpp::wrap(dmixlike_cpp(x_mat, s_mat, v_mat, pi_vec, return_log));
    return rcpp_result_gen;
END_RCPP
}
