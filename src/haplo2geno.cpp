#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "ms_files.h"
#include "gaston/matrix4.h"

using namespace Rcpp;

//[[Rcpp::export]]
XPtr<matrix4> haplo2geno(XPtr<matrix4> pA, IntegerVector H1, IntegerVector H2) {

  int nb_snps = pA->nrow;
  int nb_inds = H1.size();
  if(H2.size() != nb_inds) stop("Dimensions mismatch\n");

  XPtr<matrix4> pX(new matrix4(nb_snps, nb_inds));  

  for(int i = 0; i < nb_snps; i++) {
    for(int j = 0; j < nb_inds; j++) {
      int hap1 = H1[j] - 1;
      int hap2 = H2[j] - 1;
      char g1 = pA->get(i, hap1);
      char g2 = pA->get(i, hap2);
      if(g1 > 1 || g2 > 1) // NAs (ou autre chose que 0, 1)
        pX->set(i, j, 3);
      else
        pX->set(i, j, g1+g2);
     }
   }

  return pX;
}

RcppExport SEXP haplo2geno(SEXP pASEXP, SEXP H1SEXP, SEXP H2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type H1(H1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type H2(H2SEXP);
    rcpp_result_gen = Rcpp::wrap(haplo2geno(pA, H1, H2));
    return rcpp_result_gen;
END_RCPP
}


