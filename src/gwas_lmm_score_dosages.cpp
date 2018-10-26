#include <Rcpp.h>
#include "snp_filler_dosages.h"
#include "gwas_lmm_score.h"

using namespace Rcpp;

// [[Rcpp::export]]
List GWAS_lmm_score_dosages(CharacterVector filename, NumericVector PY, NumericMatrix P, int beg, int end) {
  snp_filler_dosages<double> S(filename, beg, end, PY.size());
  gwas_lmm_score<double> x(PY, P, S);
  x.run_tests();

  List R;
  R["id"] = wrap(S.SNP_ID);
  R["chr"] = wrap(S.CHR);
  R["pos"] = wrap(S.POS);
  R["A1"] = wrap(S.AL1);
  R["A2"] = wrap(S.AL2);
  R["freq1"] = wrap(S.F1);
  R["freq2"] = wrap(S.F2);
  R["score"] = S.L["score"];
  
  return R;
}

RcppExport SEXP gg_GWAS_lmm_score_dosages(SEXP filenameSEXP, SEXP PYSEXP, SEXP PSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_lmm_score_dosages(filename, PY, P, beg, end));
    return rcpp_result_gen;
END_RCPP
}

