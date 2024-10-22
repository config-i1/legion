// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// vSimulatorWrap
RcppExport SEXP vSimulatorWrap(SEXP arrayStates, SEXP arrayErrors, SEXP arrayF, SEXP arrayW, SEXP arrayG, SEXP modelLags);
RcppExport SEXP _legion_vSimulatorWrap(SEXP arrayStatesSEXP, SEXP arrayErrorsSEXP, SEXP arrayFSEXP, SEXP arrayWSEXP, SEXP arrayGSEXP, SEXP modelLagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type arrayStates(arrayStatesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type arrayErrors(arrayErrorsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type arrayF(arrayFSEXP);
    Rcpp::traits::input_parameter< SEXP >::type arrayW(arrayWSEXP);
    Rcpp::traits::input_parameter< SEXP >::type arrayG(arrayGSEXP);
    Rcpp::traits::input_parameter< SEXP >::type modelLags(modelLagsSEXP);
    rcpp_result_gen = Rcpp::wrap(vSimulatorWrap(arrayStates, arrayErrors, arrayF, arrayW, arrayG, modelLags));
    return rcpp_result_gen;
END_RCPP
}
// vFitterWrap
RcppExport SEXP vFitterWrap(arma::mat const& matrixY, arma::mat matrixV, arma::sp_mat& matrixF, arma::sp_mat& matrixW, arma::sp_mat& matrixG, arma::uvec& lags, char const& E, char const& T, char const& S, arma::sp_mat& matrixO);
RcppExport SEXP _legion_vFitterWrap(SEXP matrixYSEXP, SEXP matrixVSEXP, SEXP matrixFSEXP, SEXP matrixWSEXP, SEXP matrixGSEXP, SEXP lagsSEXP, SEXP ESEXP, SEXP TSEXP, SEXP SSEXP, SEXP matrixOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type matrixY(matrixYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type matrixV(matrixVSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type matrixF(matrixFSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type matrixW(matrixWSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type matrixG(matrixGSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< char const& >::type E(ESEXP);
    Rcpp::traits::input_parameter< char const& >::type T(TSEXP);
    Rcpp::traits::input_parameter< char const& >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type matrixO(matrixOSEXP);
    rcpp_result_gen = Rcpp::wrap(vFitterWrap(matrixY, matrixV, matrixF, matrixW, matrixG, lags, E, T, S, matrixO));
    return rcpp_result_gen;
END_RCPP
}
// vForecasterWrap
RcppExport SEXP vForecasterWrap(arma::mat matrixV, arma::sp_mat const& matrixF, arma::sp_mat const& matrixW, unsigned int const& nSeries, unsigned int const& hor, char const& E, char const& T, char const& S, arma::uvec& lags);
RcppExport SEXP _legion_vForecasterWrap(SEXP matrixVSEXP, SEXP matrixFSEXP, SEXP matrixWSEXP, SEXP nSeriesSEXP, SEXP horSEXP, SEXP ESEXP, SEXP TSEXP, SEXP SSEXP, SEXP lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type matrixV(matrixVSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat const& >::type matrixF(matrixFSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat const& >::type matrixW(matrixWSEXP);
    Rcpp::traits::input_parameter< unsigned int const& >::type nSeries(nSeriesSEXP);
    Rcpp::traits::input_parameter< unsigned int const& >::type hor(horSEXP);
    Rcpp::traits::input_parameter< char const& >::type E(ESEXP);
    Rcpp::traits::input_parameter< char const& >::type T(TSEXP);
    Rcpp::traits::input_parameter< char const& >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type lags(lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(vForecasterWrap(matrixV, matrixF, matrixW, nSeries, hor, E, T, S, lags));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_legion_vSimulatorWrap", (DL_FUNC) &_legion_vSimulatorWrap, 6},
    {"_legion_vFitterWrap", (DL_FUNC) &_legion_vFitterWrap, 10},
    {"_legion_vForecasterWrap", (DL_FUNC) &_legion_vForecasterWrap, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_legion(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
