#include <Rcpp.h>
#include "matrix-varia.h"
#include "matrix4.h"

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

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

template<typename scalar_t>
List GWAS_lmm_score(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end) {
  Rcout << "coucou\n";
  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");


  // copy has to be done when sclar_t = float
  VECTOR<scalar_t> Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  MATRIX<scalar_t> PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
  // ***


  int r = end-beg+1;

  VECTOR<scalar_t> SNP(n);
  NumericVector s(r);
  double t, v;  
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      s(i-beg) = NAN;
      continue;
    }

    fill_SNP<scalar_t>(pA->data[i], pA->true_ncol, pA->ncol, SNP, (scalar_t) p(i));

    //    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP); // marche pas
    //    v = (PP*SNP).dot(SNP);                          // marche mais n'utilise pas la symétrie
    v = (PP.template selfadjointView<Lower>() * SNP).dot(SNP);
    
    t = SNP.dot(Py);
    s(i-beg) = t*t/v;
  }
  
  List S;
  S["score"] = s;

  return S;
}


template<typename scalar_t>
List dominant_GWAS_lmm_score(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end) {
  Rcout << "coucou\n";
  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");


  // copy has to be done when sclar_t = float
  VECTOR<scalar_t> Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  MATRIX<scalar_t> PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
  // ***


  int r = end-beg+1;

  VECTOR<scalar_t> SNP(n);
  NumericVector s(r);
  double t, v;  
  
  for(int i = beg; i <= end; i++) {
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      s(i-beg) = NAN;
      continue;
    }

    fill_SNP_dominant<scalar_t>(pA->data[i], pA->true_ncol, pA->ncol, SNP, (scalar_t) p(i));

    //    v = (PP.selfadjointView<Lower>()*SNP).dot(SNP); // marche pas
    //    v = (PP*SNP).dot(SNP);                          // marche mais n'utilise pas la symétrie
    v = (PP.template selfadjointView<Lower>() * SNP).dot(SNP);
    
    t = SNP.dot(Py);
    s(i-beg) = t*t/v;
  }
  
  List S;
  S["score"] = s;

  return S;
}


RcppExport SEXP gg_GWAS_lmm_score(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP );
	Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< int >::type beg(begSEXP );
        Rcpp::traits::input_parameter< int >::type end(endSEXP );
        List __result = GWAS_lmm_score<double>(pA, PY, P, mu, beg, end);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_GWAS_lmm_score_f(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_lmm_score<float>(pA, PY, P, mu, beg, end));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_dominant_GWAS_lmm_score(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(dominant_GWAS_lmm_score<double>(pA, PY, P, mu, beg, end));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_dominant_GWAS_lmm_score_f(SEXP pASEXP, SEXP PYSEXP, SEXP PSEXP, SEXP muSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(dominant_GWAS_lmm_score<float>(pA, PY, P, mu, beg, end));
    return rcpp_result_gen;
END_RCPP
}

