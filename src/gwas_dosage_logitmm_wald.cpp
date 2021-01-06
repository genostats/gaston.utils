#include <Rcpp.h>
#include <string>
#include "gaston/matrix-varia.h"
#include "gaston/ai-reml-logit-1k-covar.h"
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"
#include <cmath>

//[[Rcpp::export]]
List GWAS_dosage_logitmm_wald_f(CharacterVector filename, NumericVector Y, NumericMatrix X, NumericMatrix K, int beg, int end, double tol) {

  dosages in(filename);
  int n = Y.size();
  int r = X.ncol();

  // recopiage des matrices... en float
  MatrixXf y(n,1);
  MatrixXf x(n,r);
  MatrixXf kk(n,n);
  for(int i = 0; i < n; i++) y(i,0) = (float) Y[i];

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < r; j++)
      x(i,j) = (float) X(i,j);

  for(int i = 0; i < n; i++) 
    for(int j = 0; j < n; j++)
      kk(i,j) = (float) K(i,j);

  // declare vectors containing result
  VectorXf TAU(end-beg+1);
  VectorXf BETA(end-beg+1);
  VectorXf SDBETA(end-beg+1);

  // initial values for beta, tau
  float tau = 0; 
  int niter;
  MatrixXf P(n,n);
  VectorXf omega(n);
  VectorXf beta(r);
  MatrixXf varbeta(r,r);

  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> SNP_ID, CHR, AL1, AL2;
  std::vector<int> POS;
  std::vector<float> dosage;
  std::vector<double> F1;
  std::vector<double> F2;

  int i = 0;
  while( in.read_line(dosage, snp_id, snp_pos, chr, A1, A2) ) {
    i++;
    if(i < beg) {
      dosage.clear();
      continue;
    }
    if(i > end) {
      dosage.clear();
      break;
    }
    SNP_ID.push_back(snp_id);
    POS.push_back(snp_pos);
    CHR.push_back(chr);
    AL1.push_back(A1);
    AL2.push_back(A2);

    // ajout deux cols fréquences
    double s = std::accumulate(dosage.begin(), dosage.end(), 0.0)/dosage.size()/2.0;
    F1.push_back(1.0 - s);
    F2.push_back(s);

    // remplir dernière colonne de x par dosage (centré)...
    for(int ii = 0; ii < n; ii++) x(ii,r-1) = dosage[ii] - 2.0*s;

    // use last computed tau as starting point...
    if( std::isnan(tau) ) tau = 0;
    AIREML1_logit_f(y, x, kk, true, 0, 50, tol, false, tau, niter, P, omega, beta, varbeta, true, false);

    TAU(i-beg) = tau;
    BETA(i-beg) = beta(r-1);
    SDBETA(i-beg) = sqrt(varbeta(r-1,r-1));
    
    dosage.clear();
  }

  List R;
  R["id"] = wrap(SNP_ID);
  R["chr"] = wrap(CHR);
  R["pos"] = wrap(POS);
  R["A1"] = wrap(AL1);
  R["A2"] = wrap(AL2);
  R["freq1"] = wrap(F1);
  R["freq2"] = wrap(F2);
  R["tau"] = TAU;
  R["beta"] = BETA;
  R["sd"] = SDBETA;
  return R;
}

RcppExport SEXP GWAS_dosage_logitmm_wald_f(SEXP filenameSEXP, SEXP YSEXP, SEXP XSEXP, SEXP KSEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_dosage_logitmm_wald_f(filename, Y, X, K, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

