#include <RcppArmadillo.h>
#include "gaston/matrix4.h"
#include "gaston/logit.h"
#include <iostream>

using namespace Rcpp;

// NOTE : cette version est bcp plus rapide. Pourquoi ? C'est un peu mystérieux, peut-être l'inversion de matrice ?!
void my_logistic_model(const arma::mat & y, const arma::mat & x, double eps, arma::vec & beta, arma::vec & W, arma::mat & XWX_i) {
  int n(y.n_rows), p(x.n_cols);
  arma::mat XWX(p,p), WX(n,p);
  arma::vec pi(n);
  arma::vec U(p);
  double d, log_d;
  double U_norm(1);
 
  W.set_size(n); // nécessaire parce que ... 
  W.zeros();
  beta.zeros();
  while (U_norm > eps) {
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - dot(x.row(j), beta) ) );
      W(j) = pi(j)*(1-pi(j)); 
    }
    U = x.t() *(y - pi);

    WX = diagmat(W) * x;
    XWX = x.t() * WX;

    // sym_inverse(XWX, XWX_i, log_d, d, 1e-5);
    // if(std::abs(d) < 1e-5) {
    
    if(!inv_sympd(XWX_i, XWX)) { // matrix is singular
      for(int i = 0; i < p; i++) {
        beta(i) = NAN;
        for(int j = 0; j < p; j++) XWX_i(i,j) = NAN;
      }
      return;
    }

    beta += XWX_i*U;
    U_norm = norm(U);
  }
}

//[[Rcpp::export]]
NumericMatrix fit_allelic_freq(XPtr<matrix4> pA, NumericVector p, arma::mat & x,
                       int beg, int end, double tol) {
  int r = x.n_cols;
  int n = x.n_rows;
  if(n != pA->ncol)
    stop("Dimensions mismatch");

  // for the result
  NumericMatrix BETA(r,end-beg+1);

  arma::vec beta(r), y(n);
  arma::vec W;
  arma::mat varbeta(r,r);
  beta.zeros();
  for(int i = beg; i <= end; i++) {
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      BETA(i-beg,_) = rep(NAN, r);
      continue;
     }
    // remplir y [le hack pour les données manquantes est pas terrible, on ferait mieux de virer l'individu, mais bon... ça complique]
    std::vector<unsigned int> indices;

    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        if( (xx&3) != 3 ) {
          y(4*ii+ss) = (0.5 * (double) (xx&3));
          indices.push_back(4*ii+ss);
        }
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        if( (xx&3) != 3 ) {
          y(4*ii+ss) = (0.5 * (double) (xx&3));
          indices.push_back(4*ii+ss);
        }
        xx >>= 2;
      }
    }
    // ------
    arma::uvec uvec_indices( &indices[0], indices.size(), false, true); // vector uvec sans copie de la mémoire
    my_logistic_model(y.rows(uvec_indices), x.rows(uvec_indices), tol, beta, W, varbeta);
    for(int j = 0; j < r; j++)
      BETA(j,i-beg) = beta(j);
  }
  return BETA;
}

//[[Rcpp::export]]
NumericVector fst(XPtr<matrix4> pA, NumericVector p, arma::mat & x,
                  int beg, int end, double tol) {
  int r = x.n_cols;
  int n = x.n_rows;
  if(n != pA->ncol)
    stop("Dimensions mismatch");

  // for the result
  NumericVector FST(end-beg+1);

  arma::vec beta(r), y(n);
  arma::vec W;
  arma::mat varbeta(r,r);
  beta.zeros();
  for(int i = beg; i <= end; i++) {
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      FST(i-beg) = NAN;
      continue;
     }
    // remplir y [le hack pour les données manquantes est pas terrible, on ferait mieux de virer l'individu, mais bon... ça complique]
    std::vector<unsigned int> indices;

    for(int ii = 0; ii < pA->true_ncol-1; ii++) {
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4; ss++) {
        if( (xx&3) != 3 ) {
          y(4*ii+ss) = (0.5 * (double) (xx&3));
          indices.push_back(4*ii+ss);
        }
        xx >>= 2;
      }
    }
    { int ii = pA->true_ncol-1;
      uint8_t xx = pA->data[i][ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < pA->ncol; ss++) {
        if( (xx&3) != 3 ) {
          y(4*ii+ss) = (0.5 * (double) (xx&3));
          indices.push_back(4*ii+ss);
        }
        xx >>= 2;
      }
    }
    // ------
    arma::uvec uvec_indices( &indices[0], indices.size(), false, true); // vector uvec sans copie de la mémoire
    my_logistic_model(y.rows(uvec_indices), x.rows(uvec_indices), tol, beta, W, varbeta);
    if(i == 0) Rcout << W.t() << "\n";
    FST(i-beg) = 1.0 - mean(W)/(p(i)*(1-p(i)));
  }
  return FST;
}   


RcppExport SEXP gg_fit_allelic_freq(SEXP pASEXP, SEXP pSEXP, SEXP xSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_allelic_freq(pA, p, x, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_fst(SEXP pASEXP, SEXP muSEXP, SEXP XSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(XSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(fst(pA, mu, x, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}
