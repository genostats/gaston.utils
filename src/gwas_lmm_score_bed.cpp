#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "gwas_lmm_score.h"

List GWAS_lmm_score_bed(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end) {
  snp_filler_additive_bed<double> S(pA, p, beg, end);
  gwas_lmm_score<double> x(PY, P, S);
  x.run_tests();
  return S.L;
}

List dominant_GWAS_lmm_score_bed(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end) {
  snp_filler_dominant_bed<double> S(pA, p, beg, end);
  gwas_lmm_score<double> x(PY, P, S);
  x.run_tests();
  return S.L;
}


RcppExport SEXP gg_GWAS_lmm_score_bed(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP pSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_lmm_score_bed(pA, PY, P, p, beg, end));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_dominant_GWAS_lmm_score_bed(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP pSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(dominant_GWAS_lmm_score_bed(pA, PY, P, p, beg, end));
    return rcpp_result_gen;
END_RCPP
}


