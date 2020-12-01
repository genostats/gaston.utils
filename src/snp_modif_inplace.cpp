#include <Rcpp.h>
#include "gaston/matrix4.h"
using namespace Rcpp;

inline void set_snp(uint8_t * data, size_t j, uint8_t val) {
  uint8_t & a = data[j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  a |= (val << ((j%4)*2)); // set to val
}

//[[Rcpp::export]]
void snp_modif_in_place(XPtr<matrix4> p_A, size_t snp, NumericVector val) {
  if(snp >= p_A->nrow) stop("SNP index out of range");
  if(val.length() != p_A->ncol) stop("size mismatch");
  uint8_t * d = p_A->data[snp];
  for(size_t j = 0 ; j < p_A->ncol; j++) {
   set_snp(d, j, val[j]);
  }
}

RcppExport SEXP snp_modif_in_place(SEXP p_ASEXP, SEXP snpSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< size_t >::type snp(snpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type val(valSEXP);
    snp_modif_in_place(p_A, snp, val);
    return R_NilValue;
END_RCPP
}

