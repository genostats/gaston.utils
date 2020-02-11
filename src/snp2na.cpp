#include <Rcpp.h>
#include "gaston/matrix4.h"
using namespace Rcpp;

void set2na(uint8_t * data, size_t j) {
  uint8_t & a = data[j/4];
  a |= (3 << ((j%4)*2)); // set to 3 = NA
}

//[[Rcpp::export]]
void snp2na(XPtr<matrix4> p_A, size_t snp, LogicalVector w) {
  if(snp >= p_A->nrow) stop("SNP index out of range");
  uint8_t * d = p_A->data[snp];
  for(size_t j = 0 ; j < p_A->ncol; j++) {
    if(w[j]) set2na(d, j);
  }
}

RcppExport SEXP snp2na(SEXP p_ASEXP, SEXP snpSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< size_t >::type snp(snpSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type w(wSEXP);
    snp2na(p_A, snp, w);
    return R_NilValue;
END_RCPP
}


