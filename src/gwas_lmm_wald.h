#include <Rcpp.h>
#include <RcppEigen.h>
#include "diago3.h"
#include "snp_filler.h"
#include <cmath>

#ifndef GWAS_LMM_WALD
#define GWAS_LMM_WALD
template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class gwas_lmm_wald {
  public:
  int n, r, p;
  scalar_t tol;
  VECTOR<scalar_t> sigma;
  MATRIX<scalar_t> u;
  snp_filler<scalar_t> & S;

  MATRIX<scalar_t> x;
  VECTOR<scalar_t> y;

  gwas_lmm_wald(NumericVector Y, NumericMatrix X,
                int p_, NumericVector Sigma, NumericMatrix U, double tol_, snp_filler<scalar_t> & S_)
  : n(Sigma.size()), r(X.ncol()), p(p_), tol(tol_), sigma(n), u(n,n), S(S_) {

    if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n) 
      stop("Dimensions mismatch");

    // copie sigma et u [ en double on pourrait faire des map matrix... ]
    for(int i = 0; i < n; i++) 
      sigma[i] = Sigma[i];

    for(int j = 0; j < n; j++) 
      for(int i = 0; i < n; i++)
        u(i,j) = U(i,j);

     // copie matrices X et Y avant produit (nécessaire en float)
    VECTOR<scalar_t> y0(n);
    for(int i = 0; i < n; i++)
      y0[i] = Y[i];

    MATRIX<scalar_t> x0(n, r);
    for(int j = 0; j < r; j++)
      for(int i = 0; i < n; i++)
        x0(i,j) = X(i,j);

   // --------------------------------------------
    x = u.transpose() * x0;
    y = u.transpose() * y0;
  }

  void run_tests() {

    VECTOR<scalar_t> SNP(n);

    // declare vectors containing result
    std::vector<double> H2, BETA, SDBETA;

    // object for likelihood maximization
    diag_lmm_likelihood<scalar_t> A(p, y, x, sigma);

    scalar_t h2 = 0;

    while( S.snp_fill(&SNP[0]) ) {

      A.X.col(r-1) = u.transpose() * SNP;

      // likelihood maximization
      h2 = (h2 > 0.9)?0.9:h2;
      A.newton_max( h2, 0, 0.99, tol, 10, false); // max_iter = 10... 
    
      // CALCUL DES BLUPS 
      VECTOR<scalar_t> beta, omega;
      A.blup(h2, beta, omega, false, true);

      if(A.d != 0) {
        H2.push_back(h2);
        BETA.push_back(beta(r-1));
        SDBETA.push_back(sqrt(A.v*A.XViX_i(r-1,r-1)));
      } else {
        H2.push_back(NAN);
        BETA.push_back(NAN);
        SDBETA.push_back(NAN);
      }
    }

    S.L["h2"] = wrap(H2);
    S.L["beta"] = wrap(BETA);
    S.L["sd"] = wrap(SDBETA);
  }
};

#endif
