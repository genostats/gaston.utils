#include <Rcpp.h>
#include <RcppEigen.h>
#include "approx_pql.h"
#include "snp_filler.h"
#include <cmath>

#ifndef GWAS_LMM_WALD
#define GWAS_LMM_WALD
template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class gwas_approx_pql {
  public:
  int n, r;
  scalar_t tol;
  VECTOR<scalar_t> sigma;
  MATRIX<scalar_t> u;
  snp_filler<scalar_t> & S;

  MATRIX<scalar_t> x;
  VECTOR<scalar_t> y;

  gwas_approx_pql(NumericVector Y, NumericMatrix X,
                  NumericVector Sigma, NumericMatrix U, double tol_, snp_filler<scalar_t> & S_)
  : n(Sigma.size()), r(X.ncol()), tol(tol_), sigma(n), u(n,n), S(S_), x(n, X.ncol()), y(n) {

    if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n) 
      stop("Dimensions mismatch");

    // copie sigma et u [ en double on pourrait faire des map matrix... ]
    for(int i = 0; i < n; i++) 
      sigma[i] = Sigma[i];

    for(int j = 0; j < n; j++) 
      for(int i = 0; i < n; i++)
        u(i,j) = U(i,j);

     // copie matrices X et Y [idem]
    for(int i = 0; i < n; i++)
      y[i] = Y[i];

    for(int j = 0; j < r; j++)
      for(int i = 0; i < n; i++)
        x(i,j) = X(i,j);
  }

  void run_tests() {

    VECTOR<scalar_t> SNP(n);

    // declare vectors containing result
    std::vector<double> H2, BETA, SDBETA;

    // object for likelihood maximization
    approx_pql<scalar_t> A(y, x, sigma, u, tol);

    scalar_t h2 = 0;
    while( S.snp_fill(&A.X(0,r-1)) ) {

      // optimisation
      h2 = (h2 > 0.9)?0.9:h2;
      A.optimize( h2, false); 
    
  //  if(A.dlmm_object.d != 0) {
        H2.push_back(A.h2);
        BETA.push_back(A.Beta(r-1));
        SDBETA.push_back(sqrt(A.VarBeta(r-1,r-1)));
  //  } else {
  //    H2.push_back(NAN);
  //    BETA.push_back(NAN);
  //    SDBETA.push_back(NAN);
  //  }
    }

    S.L["h2"] = wrap(H2);
    S.L["beta"] = wrap(BETA);
    S.L["sd"] = wrap(SDBETA);
  }
};

#endif
