// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// slopeOPtransfer
List slopeOPtransfer(std::vector<double> data, std::vector<double> states, double penalty, std::string constraint, double minAngle, std::string type);
RcppExport SEXP _slopeOP_slopeOPtransfer(SEXP dataSEXP, SEXP statesSEXP, SEXP penaltySEXP, SEXP constraintSEXP, SEXP minAngleSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type states(statesSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< std::string >::type constraint(constraintSEXP);
    Rcpp::traits::input_parameter< double >::type minAngle(minAngleSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(slopeOPtransfer(data, states, penalty, constraint, minAngle, type));
    return rcpp_result_gen;
END_RCPP
}
// slopeSNtransfer
List slopeSNtransfer(std::vector<double> data, std::vector<double> states, unsigned int nbSegments, std::string constraint);
RcppExport SEXP _slopeOP_slopeSNtransfer(SEXP dataSEXP, SEXP statesSEXP, SEXP nbSegmentsSEXP, SEXP constraintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type states(statesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbSegments(nbSegmentsSEXP);
    Rcpp::traits::input_parameter< std::string >::type constraint(constraintSEXP);
    rcpp_result_gen = Rcpp::wrap(slopeSNtransfer(data, states, nbSegments, constraint));
    return rcpp_result_gen;
END_RCPP
}
// linearOP
List linearOP(std::vector<double> x, std::vector<double> data, double penalty, bool cc);
RcppExport SEXP _slopeOP_linearOP(SEXP xSEXP, SEXP dataSEXP, SEXP penaltySEXP, SEXP ccSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< bool >::type cc(ccSEXP);
    rcpp_result_gen = Rcpp::wrap(linearOP(x, data, penalty, cc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_slopeOP_slopeOPtransfer", (DL_FUNC) &_slopeOP_slopeOPtransfer, 6},
    {"_slopeOP_slopeSNtransfer", (DL_FUNC) &_slopeOP_slopeSNtransfer, 4},
    {"_slopeOP_linearOP", (DL_FUNC) &_slopeOP_linearOP, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_slopeOP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
