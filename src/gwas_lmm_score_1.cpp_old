#include <Rcpp.h>
#include "matrix-varia.h"
#include "matrix4.h"

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

  gwas_lmm_score(NumericVector PY, NumericMatrix P) 
  : n(PY.size()), Py(n), PP(n,n), SNP(n) {
    if(P.nrow() != n || P.ncol() != n) 
      stop("Dimensions mismatch\n");

    // copy has to be done when scalar_t = float
    for(int i = 0; i < n; i++) Py(i) = PY[i];

    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) 
        PP(i,j) = P(i,j);
  }

  virtual bool snp_fill() = 0;
  virtual bool current_snp_monomorphic() = 0;

  List run_tests() {
    std::vector<double> s;
    double t, v;  
    while( snp_fill() ) {
      if( current_snp_monomorphic() ) {
        s.push_back(NAN);
        continue;
      }
      //    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP); // marche pas
      //    v = (PP*SNP).dot(SNP);                          // marche mais n'utilise pas la symétrie
      v = (PP.template selfadjointView<Lower>() * SNP).dot(SNP); // la solution !
    
      t = SNP.dot(Py);
      s.push_back(t*t/v);
    }
  
    List S;
    S["score"] = wrap(s);

    return S;
  }
};

template<typename scalar_t>
class gwas_lmm_score_additive_bed : public gwas_lmm_score<scalar_t> {
  public:
  XPtr<matrix4> pA;  
  int ncol, true_ncol;
  NumericVector p; 
  int beg, end;
  int i;
  bool monomorphic;
  gwas_lmm_score_additive_bed(NumericVector PY, NumericMatrix P, XPtr<matrix4> pA_, NumericVector p_, int beg_, int end_)
    : gwas_lmm_score<scalar_t>(PY, P), pA(pA_), ncol(pA->ncol), true_ncol(pA->true_ncol), p(p_), beg(beg_), end(end_), i(beg), monomorphic(true) {};

  bool current_snp_monomorphic() {
    return monomorphic;
  }

  bool snp_fill() {
    if(i > end) {
      monomorphic = true; 
      return false; 
    }
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      monomorphic = true;
      i++;
      return true;
    }
    uint8_t * snp = pA-> data[i];
    scalar_t mu = 2*p(i);
    for(int ii = 0; ii < true_ncol-1; ii++) {
      uint8_t x = snp[ii];
      for(int ss = 0; ss < 4; ss++) {
        this->SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu; // le 'this' est nécessaire pour que le compilateur sache d'où vient SNP
        x >>= 2;
      }
    }
    { int ii = true_ncol-1;
      uint8_t x = snp[ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
        this->SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu;
        x >>= 2;
      }
    }
    i++;
    monomorphic = false;
    return true;
  }
};

//[[Rcpp::export]]
List GWAS_lmm_score_1(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end) {
  gwas_lmm_score_additive_bed<double> x(PY, P, pA, p, beg, end);
  return x.run_tests();
}


/*
template<typename scalar_t>
void fill_SNP(uint8_t * snp, int true_ncol, int ncol, VECTOR<scalar_t> & SNP, scalar_t p) {
  // récupérer SNP
  scalar_t mu = 2*p;
  for(int ii = 0; ii < true_ncol-1; ii++) {
    uint8_t x = snp[ii];
    for(int ss = 0; ss < 4; ss++) {
      SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu;
      x >>= 2;
    }
  }
  { int ii = true_ncol-1;
    uint8_t x = snp[ii];
    for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
      SNP(4*ii+ss) = ((x&3) != 3)?(x&3):mu;
      x >>= 2;
    }
  }
}

template<typename scalar_t>
void fill_SNP_dominant(uint8_t * snp, int true_ncol, int ncol, VECTOR<scalar_t> & SNP, scalar_t p) {
  // récupérer SNP
  scalar_t h[4];
  h[0] = p / (1-p);
  h[1] = -1;
  h[2] = (1-p) / p;
  h[3] = 0; // misson values imputed as 0
  for(int ii = 0; ii < true_ncol-1; ii++) {
    uint8_t x = snp[ii];
    for(int ss = 0; ss < 4; ss++) {
      SNP(4*ii+ss) = h[x&3];
      x >>= 2;
    }
  }
  { int ii = true_ncol-1;
    uint8_t x = snp[ii];
    for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
      SNP(4*ii+ss) = h[x&3];
      x >>= 2;
    }
  }
}
*/

RcppExport SEXP gg_GWAS_lmm_score_1(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP pSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_lmm_score_1(pA, PY, P, p, beg, end));
    return rcpp_result_gen;
END_RCPP
}


