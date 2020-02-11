#include <Rcpp.h>
#include "snp_filler_dosages.h"
#include "gwas_logit_wald.h"

using namespace Rcpp;

// [[Rcpp::export]]
List GWAS_logit_wald_dosages(CharacterVector filename, NumericVector Y, NumericMatrix X, int beg, int end, double tol) {
  snp_filler_dosages<float> S(filename, beg, end, Y.size());
  gwas_logit_wald<float> x(Y, X, tol, S);
  x.run_tests();

  List R;
  R["id"] = wrap(S.SNP_ID);
  R["chr"] = wrap(S.CHR);
  R["pos"] = wrap(S.POS);
  R["A1"] = wrap(S.AL1);
  R["A2"] = wrap(S.AL2);
  R["freq1"] = wrap(S.F1);
  R["freq2"] = wrap(S.F2);
  R["beta"] = S.L["beta"];
  R["sd"] = S.L["sd"];
  
  return R;
}



RcppExport SEXP gg_GWAS_logit_wald_dosages(SEXP filenameSEXP, SEXP YSEXP, SEXP XSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_logit_wald_dosages(filename, Y, X, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

