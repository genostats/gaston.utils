#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include "logit_model.h"
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"
#include <cmath>

#ifndef GWASLOGITWALD
#define GWASLOGITWALD
template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class gwas_logit_wald {
  public:
  int n, r;
  scalar_t tol;
  VECTOR<scalar_t> y;
  MATRIX<scalar_t> x;
  snp_filler<scalar_t> & S;

  gwas_logit_wald(NumericVector Y, NumericMatrix X, double tol_, snp_filler<scalar_t> & S_) : 
    n(Y.size()), r(X.ncol()), tol(tol_), y(n), x(n,r), S(S_) {
    // recopiage des matrices... nécessaire en float
    for(int i = 0; i < n; i++) y(i) = (scalar_t) Y[i];

    for(int i = 0; i < n; i++) 
      for(int j = 0; j < r; j++)
        x(i,j) = (scalar_t) X(i,j);

  }

  void run_tests() {
    std::vector<double> BETA, SDBETA;

    VECTOR<scalar_t> beta(r); 
    beta.setZero();

    // remplir dernière colonne de x 
    MATRIX<scalar_t>  varbeta(r,r);
    while( S.snp_fill( &x(00,r-1) )) {
      if( S.current_snp_monomorphic() ) { // pas la peine d'aller plus loin
        BETA.push_back(NAN);
        SDBETA.push_back(NAN);
        continue;
      }
      // verifier si la matrice est singuliere
      scalar_t d = (x.transpose() * x).determinant();
      if( d < 1e-4 ) { // seuil très arbitraire (devrait être ok si on a fait la dec QR de X avant)
        BETA.push_back(NAN);
        SDBETA.push_back(NAN);
      } else {
        logistic_model2<scalar_t>(y, x, beta, varbeta, tol);
        BETA.push_back(beta(r-1));
        SDBETA.push_back(sqrt(varbeta(r-1,r-1)));
      }
    }

    S.L["beta"] = wrap(BETA);
    S.L["sd"] = wrap(SDBETA);
  }
};


#endif
