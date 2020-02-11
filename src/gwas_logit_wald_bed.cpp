#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "gwas_logit_wald.h"

//[[Rcpp::export]]
List GWAS_logit_wald_bed(XPtr<matrix4> pA, NumericVector p, NumericVector Y, NumericMatrix X,
                       int beg, int end, double tol) {
  snp_filler_additive_bed<float> S(pA, p, beg, end);
  gwas_logit_wald<float> x(Y, X, tol, S);
  x.run_tests();
  return S.L;
}

RcppExport SEXP gg_GWAS_logit_wald_bed(SEXP pASEXP, SEXP pSEXP, SEXP YSEXP, SEXP XSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_logit_wald_bed(pA, p, Y, X, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

