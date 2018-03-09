#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"

using namespace Rcpp;
using namespace RcppParallel;


struct para_w_colmean : public Worker {
  const matrix4 & A;
  const RVector<double> alpha;
  const size_t ncol;
  const size_t true_ncol;
  const double * v;
  
  //output
  double * means;

  //constructeur
  para_w_colmean(matrix4 & A_, NumericVector alpha_, double * means_)
           : A(A_), alpha(alpha_), ncol(A.ncol), true_ncol(A.true_ncol), means(means_) {}

  //worker
  void operator()(size_t beg, size_t end) {
    for(size_t i = beg; i < end; i++) {
      double S = 0, w = 0;
      size_t k = 0;
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t x = A.data[i][j];
        for(int ss = 0; ss < 4 && 4*j + ss < ncol; ss++) {
          if((x&3) != 3) {
            S += alpha[k]*((double) (x&3));
            w += alpha[k++];
          }
          x >>= 2;
        }
      }
      means[i] = S/w;
    }
  }
};

// [[Rcpp::export]]
NumericVector weighted_col_means(XPtr<matrix4> p_A, NumericVector alpha) {
  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds
  if(m != alpha.size()) stop("Dimensions mismatch");

  NumericVector R(n);
  para_w_colmean X(*p_A, alpha, R.begin());
  parallelFor(0, n, X, 100);
  return R;
}

RcppExport SEXP gg_weighted_col_means(SEXP p_ASEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_col_means(p_A, alpha));
    return rcpp_result_gen;
END_RCPP
}

