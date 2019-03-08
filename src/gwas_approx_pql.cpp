#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "gwas_approx_pql.h"

//[[Rcpp::export]]
List GWAS_approx_pql_bed(XPtr<matrix4> pA, NumericVector pp, NumericVector Y, NumericMatrix X,
                       NumericVector Sigma, NumericMatrix U, int beg, int end, double tol) {
  snp_filler_additive_bed<double> S(pA, pp, beg, end);
  gwas_approx_pql<double> x(Y, X, Sigma, U, tol, S);
  x.run_tests();
  return S.L;
}

RcppExport SEXP gg_GWAS_approx_pql_bed(SEXP pASEXP, SEXP ppSEXP, SEXP YSEXP, SEXP XSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_approx_pql_bed(pA, pp, Y, X, Sigma, U, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

