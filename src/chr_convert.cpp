#include <Rcpp.h>
#include "chr_convert.h"

using namespace Rcpp;

void set_chr_ids(List L) {
  chr_ids = L;
}

RcppExport SEXP gg_set_chr_ids(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    set_chr_ids(L);
    return R_NilValue;
END_RCPP
}



