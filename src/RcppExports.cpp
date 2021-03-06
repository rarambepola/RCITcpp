// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// cov_cpp2
Eigen::MatrixXd cov_cpp2(Eigen::MatrixXd x, Eigen::MatrixXd y);
RcppExport SEXP _RCITcpp_cov_cpp2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cov_cpp2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// normalise_cpp2
Eigen::MatrixXd normalise_cpp2(Eigen::MatrixXd x);
RcppExport SEXP _RCITcpp_normalise_cpp2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(normalise_cpp2(x));
    return rcpp_result_gen;
END_RCPP
}
// make_rff_cpp
Eigen::MatrixXd make_rff_cpp(Eigen::MatrixXd x, Eigen::MatrixXd w, double sigma, Eigen::VectorXd b);
RcppExport SEXP _RCITcpp_make_rff_cpp(SEXP xSEXP, SEXP wSEXP, SEXP sigmaSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(make_rff_cpp(x, w, sigma, b));
    return rcpp_result_gen;
END_RCPP
}
// RIT_cpp
double RIT_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y, double median_x, double median_y, Eigen::MatrixXd w, Eigen::VectorXd b, bool return_ts, double n_bs);
RcppExport SEXP _RCITcpp_RIT_cpp(SEXP xSEXP, SEXP ySEXP, SEXP median_xSEXP, SEXP median_ySEXP, SEXP wSEXP, SEXP bSEXP, SEXP return_tsSEXP, SEXP n_bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type median_x(median_xSEXP);
    Rcpp::traits::input_parameter< double >::type median_y(median_ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ts(return_tsSEXP);
    Rcpp::traits::input_parameter< double >::type n_bs(n_bsSEXP);
    rcpp_result_gen = Rcpp::wrap(RIT_cpp(x, y, median_x, median_y, w, b, return_ts, n_bs));
    return rcpp_result_gen;
END_RCPP
}
// RCIT_cpp
double RCIT_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd z, double median_x, double median_y, double median_z, Eigen::MatrixXd w, Eigen::MatrixXd wy, Eigen::MatrixXd wz, Eigen::VectorXd b, Eigen::VectorXd by, Eigen::VectorXd bz, bool return_ts, double n_bs);
RcppExport SEXP _RCITcpp_RCIT_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP median_xSEXP, SEXP median_ySEXP, SEXP median_zSEXP, SEXP wSEXP, SEXP wySEXP, SEXP wzSEXP, SEXP bSEXP, SEXP bySEXP, SEXP bzSEXP, SEXP return_tsSEXP, SEXP n_bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type median_x(median_xSEXP);
    Rcpp::traits::input_parameter< double >::type median_y(median_ySEXP);
    Rcpp::traits::input_parameter< double >::type median_z(median_zSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type wy(wySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type wz(wzSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type by(bySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type bz(bzSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ts(return_tsSEXP);
    Rcpp::traits::input_parameter< double >::type n_bs(n_bsSEXP);
    rcpp_result_gen = Rcpp::wrap(RCIT_cpp(x, y, z, median_x, median_y, median_z, w, wy, wz, b, by, bz, return_ts, n_bs));
    return rcpp_result_gen;
END_RCPP
}
// RCIT_disag_cpp
double RCIT_disag_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd z, double median_x, double median_y, double median_z, Eigen::MatrixXd w, Eigen::MatrixXd wy, Eigen::MatrixXd wz, Eigen::VectorXd b, Eigen::VectorXd by, Eigen::VectorXd bz, IntegerVector polygon_start_index, IntegerVector polygon_end_index, bool population_supplied, NumericVector population, bool return_ts, double n_bs);
RcppExport SEXP _RCITcpp_RCIT_disag_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP median_xSEXP, SEXP median_ySEXP, SEXP median_zSEXP, SEXP wSEXP, SEXP wySEXP, SEXP wzSEXP, SEXP bSEXP, SEXP bySEXP, SEXP bzSEXP, SEXP polygon_start_indexSEXP, SEXP polygon_end_indexSEXP, SEXP population_suppliedSEXP, SEXP populationSEXP, SEXP return_tsSEXP, SEXP n_bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type median_x(median_xSEXP);
    Rcpp::traits::input_parameter< double >::type median_y(median_ySEXP);
    Rcpp::traits::input_parameter< double >::type median_z(median_zSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type wy(wySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type wz(wzSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type by(bySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type bz(bzSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_start_index(polygon_start_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_end_index(polygon_end_indexSEXP);
    Rcpp::traits::input_parameter< bool >::type population_supplied(population_suppliedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ts(return_tsSEXP);
    Rcpp::traits::input_parameter< double >::type n_bs(n_bsSEXP);
    rcpp_result_gen = Rcpp::wrap(RCIT_disag_cpp(x, y, z, median_x, median_y, median_z, w, wy, wz, b, by, bz, polygon_start_index, polygon_end_index, population_supplied, population, return_ts, n_bs));
    return rcpp_result_gen;
END_RCPP
}
// RIT_disag_cpp
double RIT_disag_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y, double median_x, double median_y, Eigen::MatrixXd w, Eigen::VectorXd b, IntegerVector polygon_start_index, IntegerVector polygon_end_index, bool population_supplied, NumericVector population, bool return_ts, double n_bs);
RcppExport SEXP _RCITcpp_RIT_disag_cpp(SEXP xSEXP, SEXP ySEXP, SEXP median_xSEXP, SEXP median_ySEXP, SEXP wSEXP, SEXP bSEXP, SEXP polygon_start_indexSEXP, SEXP polygon_end_indexSEXP, SEXP population_suppliedSEXP, SEXP populationSEXP, SEXP return_tsSEXP, SEXP n_bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type median_x(median_xSEXP);
    Rcpp::traits::input_parameter< double >::type median_y(median_ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_start_index(polygon_start_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_end_index(polygon_end_indexSEXP);
    Rcpp::traits::input_parameter< bool >::type population_supplied(population_suppliedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ts(return_tsSEXP);
    Rcpp::traits::input_parameter< double >::type n_bs(n_bsSEXP);
    rcpp_result_gen = Rcpp::wrap(RIT_disag_cpp(x, y, median_x, median_y, w, b, polygon_start_index, polygon_end_index, population_supplied, population, return_ts, n_bs));
    return rcpp_result_gen;
END_RCPP
}
// RIT_disag_cpp_v2
double RIT_disag_cpp_v2(Eigen::MatrixXd x, Eigen::MatrixXd y, double median_x, double median_y, Eigen::MatrixXd w_x, Eigen::VectorXd b_x, Eigen::MatrixXd w_y, Eigen::VectorXd b_y, IntegerVector polygon_start_index, IntegerVector polygon_end_index, NumericVector population, int n_obs, bool return_ts, double n_bs);
RcppExport SEXP _RCITcpp_RIT_disag_cpp_v2(SEXP xSEXP, SEXP ySEXP, SEXP median_xSEXP, SEXP median_ySEXP, SEXP w_xSEXP, SEXP b_xSEXP, SEXP w_ySEXP, SEXP b_ySEXP, SEXP polygon_start_indexSEXP, SEXP polygon_end_indexSEXP, SEXP populationSEXP, SEXP n_obsSEXP, SEXP return_tsSEXP, SEXP n_bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type median_x(median_xSEXP);
    Rcpp::traits::input_parameter< double >::type median_y(median_ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w_x(w_xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b_x(b_xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w_y(w_ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b_y(b_ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_start_index(polygon_start_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_end_index(polygon_end_indexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ts(return_tsSEXP);
    Rcpp::traits::input_parameter< double >::type n_bs(n_bsSEXP);
    rcpp_result_gen = Rcpp::wrap(RIT_disag_cpp_v2(x, y, median_x, median_y, w_x, b_x, w_y, b_y, polygon_start_index, polygon_end_index, population, n_obs, return_ts, n_bs));
    return rcpp_result_gen;
END_RCPP
}
// RCIT_disag_cpp_v2
double RCIT_disag_cpp_v2(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd z, double median_x, double median_y, double median_z, Eigen::MatrixXd w, Eigen::MatrixXd wy, Eigen::MatrixXd wz, Eigen::VectorXd b, Eigen::VectorXd by, Eigen::VectorXd bz, IntegerVector polygon_start_index, IntegerVector polygon_end_index, NumericVector population, int n_obs, bool return_ts, double n_bs);
RcppExport SEXP _RCITcpp_RCIT_disag_cpp_v2(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP median_xSEXP, SEXP median_ySEXP, SEXP median_zSEXP, SEXP wSEXP, SEXP wySEXP, SEXP wzSEXP, SEXP bSEXP, SEXP bySEXP, SEXP bzSEXP, SEXP polygon_start_indexSEXP, SEXP polygon_end_indexSEXP, SEXP populationSEXP, SEXP n_obsSEXP, SEXP return_tsSEXP, SEXP n_bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type median_x(median_xSEXP);
    Rcpp::traits::input_parameter< double >::type median_y(median_ySEXP);
    Rcpp::traits::input_parameter< double >::type median_z(median_zSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type wy(wySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type wz(wzSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type by(bySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type bz(bzSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_start_index(polygon_start_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type polygon_end_index(polygon_end_indexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ts(return_tsSEXP);
    Rcpp::traits::input_parameter< double >::type n_bs(n_bsSEXP);
    rcpp_result_gen = Rcpp::wrap(RCIT_disag_cpp_v2(x, y, z, median_x, median_y, median_z, w, wy, wz, b, by, bz, polygon_start_index, polygon_end_index, population, n_obs, return_ts, n_bs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RCITcpp_cov_cpp2", (DL_FUNC) &_RCITcpp_cov_cpp2, 2},
    {"_RCITcpp_normalise_cpp2", (DL_FUNC) &_RCITcpp_normalise_cpp2, 1},
    {"_RCITcpp_make_rff_cpp", (DL_FUNC) &_RCITcpp_make_rff_cpp, 4},
    {"_RCITcpp_RIT_cpp", (DL_FUNC) &_RCITcpp_RIT_cpp, 8},
    {"_RCITcpp_RCIT_cpp", (DL_FUNC) &_RCITcpp_RCIT_cpp, 14},
    {"_RCITcpp_RCIT_disag_cpp", (DL_FUNC) &_RCITcpp_RCIT_disag_cpp, 18},
    {"_RCITcpp_RIT_disag_cpp", (DL_FUNC) &_RCITcpp_RIT_disag_cpp, 12},
    {"_RCITcpp_RIT_disag_cpp_v2", (DL_FUNC) &_RCITcpp_RIT_disag_cpp_v2, 14},
    {"_RCITcpp_RCIT_disag_cpp_v2", (DL_FUNC) &_RCITcpp_RCIT_disag_cpp_v2, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_RCITcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
