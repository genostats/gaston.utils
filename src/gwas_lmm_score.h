#include <Rcpp.h>
#include <RcppEigen.h>
#include "snp_filler.h"

#ifndef GWAS_LMM_SCORE
#define GWAS_LMM_SCORE

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class gwas_lmm_score {
  public:
  int n;
  VECTOR<scalar_t> Py;
  MATRIX<scalar_t> PP;
  VECTOR<scalar_t> SNP;
  snp_filler<scalar_t> & S;

  gwas_lmm_score(NumericVector PY, NumericMatrix P, snp_filler<scalar_t> & S_) 
  : n(PY.size()), Py(n), PP(n,n), SNP(n), S(S_) {
    if(P.nrow() != n || P.ncol() != n) 
      stop("Dimensions mismatch\n");

    // copy has to be done when scalar_t = float
    for(int i = 0; i < n; i++) Py(i) = PY[i];

    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) 
        PP(i,j) = P(i,j);
  }

//  virtual bool snp_fill() = 0;
//  virtual bool current_snp_monomorphic() = 0;

  void run_tests() {
    std::vector<double> s;
    double t, v;  

    while( S.snp_fill( &SNP[0] ) ) {
      if( S.current_snp_monomorphic() ) {
        s.push_back(NAN);
        continue;
      }
      //    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP); // marche pas
      //    v = (PP*SNP).dot(SNP);                          // marche mais n'utilise pas la symétrie
      v = (PP.template selfadjointView<Eigen::Lower>() * SNP).dot(SNP); // la solution !
      t = SNP.dot(Py);
      s.push_back(t*t/v);
    }

    // on ajoute le score dans S.L 
    S.L["score"] = wrap(s);
  }
};

#endif
