#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "gwas_logit_offset.h"

//[[Rcpp::export]]
List GWAS_logit_offset_bed(XPtr<matrix4> pA, NumericVector p, NumericVector Y, NumericVector Offset,
                           NumericMatrix Q, int beg, int end, double tol, int max_iter) {
  snp_filler_additive_bed<double> S(pA, p, beg, end);
  gwas_logit_offset<double> x(Y, Offset, Q, tol, max_iter, S);
  x.run_tests();
  return S.L;
}

RcppExport SEXP gg_GWAS_logit_offset_bed(SEXP pASEXP, SEXP pSEXP, SEXP YSEXP, SEXP OffsetSEXP, SEXP QSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Offset(OffsetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_logit_offset_bed(pA, p, Y, Offset, Q, beg, end, tol, max_iter));
    return rcpp_result_gen;
END_RCPP
}

