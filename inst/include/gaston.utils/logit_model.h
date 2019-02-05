#include <RcppEigen.h>
#include <cmath>
#include <iostream>
#ifndef GASTONLOGITMODEL
#define GASTONLOGITMODEL

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

using namespace Rcpp;
using namespace Eigen;

template<typename scalar_t>
void logistic_model2(const MATRIX<scalar_t> & y, const MATRIX<scalar_t> & x, scalar_t eps, VECTOR<scalar_t> & beta, MATRIX<scalar_t> & XWX_i, int max_iter) {
  int n(y.rows()), p(x.cols());
  VECTOR<scalar_t> W(n), pi(n), U(p);
  MATRIX<scalar_t> XWX(p,p), WX(n,p);

  W.setZero();

  scalar_t U_norm(1);
  beta.setZero();
  int k = 1;

  while (U_norm > eps) {
    for(int j = 0; j < n; j++) {
      pi(j) = 1/( 1 + exp( - x.row(j).dot(beta) ) );
      W(j) = pi(j)*(1-pi(j));
    }
    U.noalias() = x.transpose()*(y - pi);

    WX.noalias() = W.asDiagonal() * x;
    XWX.noalias() = x.transpose() * WX;

    beta += XWX.llt().solve(U);  // SHOULD BE BETTER THAN sym_inverse(XWX, XWX_i, log_d, d, 1e-5);

    U_norm = U.norm();
    if(k++ > max_iter) return;
  } 
  // Mais il nous faut l'inverse de XWX pour calculer l'écart type
  XWX_i = XWX.llt().solve( MatrixXd::Identity(p,p) );
}
#endif
