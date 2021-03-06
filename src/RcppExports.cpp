// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// reFit
Rcpp::List reFit(const std::vector<double>& yVec, const std::vector<double>& kVec, const int32_t& d, const int32_t& Ngen);
RcppExport SEXP _GWAlikeMeth_reFit(SEXP yVecSEXP, SEXP kVecSEXP, SEXP dSEXP, SEXP NgenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    rcpp_result_gen = Rcpp::wrap(reFit(yVec, kVec, d, Ngen));
    return rcpp_result_gen;
END_RCPP
}
// reFitR
Rcpp::List reFitR(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& repFac, const int32_t& d, const int32_t& Ngen);
RcppExport SEXP _GWAlikeMeth_reFitR(SEXP yVecSEXP, SEXP kVecSEXP, SEXP repFacSEXP, SEXP dSEXP, SEXP NgenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type repFac(repFacSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    rcpp_result_gen = Rcpp::wrap(reFitR(yVec, kVec, repFac, d, Ngen));
    return rcpp_result_gen;
END_RCPP
}
// reFitF
Rcpp::List reFitF(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<double>& xVec, const int32_t& d, const int32_t& Ngen);
RcppExport SEXP _GWAlikeMeth_reFitF(SEXP yVecSEXP, SEXP kVecSEXP, SEXP xVecSEXP, SEXP dSEXP, SEXP NgenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xVec(xVecSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    rcpp_result_gen = Rcpp::wrap(reFitF(yVec, kVec, xVec, d, Ngen));
    return rcpp_result_gen;
END_RCPP
}
// reFitRF
Rcpp::List reFitRF(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& repFac, const std::vector<double>& xVec, const int32_t& d, const int32_t& Ngen);
RcppExport SEXP _GWAlikeMeth_reFitRF(SEXP yVecSEXP, SEXP kVecSEXP, SEXP repFacSEXP, SEXP xVecSEXP, SEXP dSEXP, SEXP NgenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type repFac(repFacSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xVec(xVecSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    rcpp_result_gen = Rcpp::wrap(reFitRF(yVec, kVec, repFac, xVec, d, Ngen));
    return rcpp_result_gen;
END_RCPP
}
// gwa
Rcpp::List gwa(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwa(SEXP yVecSEXP, SEXP kVecSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwa(yVec, kVec, snps, d, Ngen, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaF
Rcpp::List gwaF(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<double>& xVec, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaF(SEXP yVecSEXP, SEXP kVecSEXP, SEXP xVecSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xVec(xVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaF(yVec, kVec, xVec, snps, d, Ngen, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaR
Rcpp::List gwaR(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& repFac, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaR(SEXP yVecSEXP, SEXP kVecSEXP, SEXP repFacSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type repFac(repFacSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaR(yVec, kVec, repFac, snps, d, Ngen, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaRF
Rcpp::List gwaRF(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& repFac, const std::vector<double> xVec, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaRF(SEXP yVecSEXP, SEXP kVecSEXP, SEXP repFacSEXP, SEXP xVecSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type repFac(repFacSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type xVec(xVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaRF(yVec, kVec, repFac, xVec, snps, d, Ngen, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaFDR
Rcpp::List gwaFDR(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nPer, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaFDR(SEXP yVecSEXP, SEXP kVecSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nPerSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nPer(nPerSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaFDR(yVec, kVec, snps, d, Ngen, nPer, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaFDRR
Rcpp::List gwaFDRR(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& repFac, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nPer, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaFDRR(SEXP yVecSEXP, SEXP kVecSEXP, SEXP repFacSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nPerSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type repFac(repFacSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nPer(nPerSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaFDRR(yVec, kVec, repFac, snps, d, Ngen, nPer, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaFDRF
Rcpp::List gwaFDRF(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<double>& xVec, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nPer, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaFDRF(SEXP yVecSEXP, SEXP kVecSEXP, SEXP xVecSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nPerSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xVec(xVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nPer(nPerSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaFDRF(yVec, kVec, xVec, snps, d, Ngen, nPer, nThr));
    return rcpp_result_gen;
END_RCPP
}
// gwaFDRRF
Rcpp::List gwaFDRRF(const std::vector<double>& yVec, const std::vector<double>& kVec, const std::vector<int32_t>& repFac, const std::vector<double>& xVec, const std::vector<int32_t>& snps, const int32_t& d, const int32_t& Ngen, const int32_t& nPer, const int32_t& nThr);
RcppExport SEXP _GWAlikeMeth_gwaFDRRF(SEXP yVecSEXP, SEXP kVecSEXP, SEXP repFacSEXP, SEXP xVecSEXP, SEXP snpsSEXP, SEXP dSEXP, SEXP NgenSEXP, SEXP nPerSEXP, SEXP nThrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type yVec(yVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type kVec(kVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type repFac(repFacSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xVec(xVecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int32_t>& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type Ngen(NgenSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nPer(nPerSEXP);
    Rcpp::traits::input_parameter< const int32_t& >::type nThr(nThrSEXP);
    rcpp_result_gen = Rcpp::wrap(gwaFDRRF(yVec, kVec, repFac, xVec, snps, d, Ngen, nPer, nThr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GWAlikeMeth_reFit", (DL_FUNC) &_GWAlikeMeth_reFit, 4},
    {"_GWAlikeMeth_reFitR", (DL_FUNC) &_GWAlikeMeth_reFitR, 5},
    {"_GWAlikeMeth_reFitF", (DL_FUNC) &_GWAlikeMeth_reFitF, 5},
    {"_GWAlikeMeth_reFitRF", (DL_FUNC) &_GWAlikeMeth_reFitRF, 6},
    {"_GWAlikeMeth_gwa", (DL_FUNC) &_GWAlikeMeth_gwa, 6},
    {"_GWAlikeMeth_gwaF", (DL_FUNC) &_GWAlikeMeth_gwaF, 7},
    {"_GWAlikeMeth_gwaR", (DL_FUNC) &_GWAlikeMeth_gwaR, 7},
    {"_GWAlikeMeth_gwaRF", (DL_FUNC) &_GWAlikeMeth_gwaRF, 8},
    {"_GWAlikeMeth_gwaFDR", (DL_FUNC) &_GWAlikeMeth_gwaFDR, 7},
    {"_GWAlikeMeth_gwaFDRR", (DL_FUNC) &_GWAlikeMeth_gwaFDRR, 8},
    {"_GWAlikeMeth_gwaFDRF", (DL_FUNC) &_GWAlikeMeth_gwaFDRF, 8},
    {"_GWAlikeMeth_gwaFDRRF", (DL_FUNC) &_GWAlikeMeth_gwaFDRRF, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_GWAlikeMeth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
