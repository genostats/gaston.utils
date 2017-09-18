#include <Rcpp.h>
#include <RcppEigen.h>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"
#include <cmath>

using namespace Rcpp;

//[[Rcpp::export]]
List GWAS_dosage_lmm_score(CharacterVector filename, NumericVector PY, NumericMatrix P, int beg, int end) {

  dosages in(filename);

  Eigen::Map<Eigen::MatrixXd> Py(as<Eigen::Map<Eigen::MatrixXd> >(PY));
  Eigen::Map<Eigen::MatrixXd> PP(as<Eigen::Map<Eigen::MatrixXd> >(P));

  int r= end-beg+1;
  int n=Py.rows();
  Eigen::VectorXd SNP(n);
  NumericVector s(r);
  double t, v;  

  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> SNP_ID, CHR, AL1, AL2;
  std::vector<int> POS;
  std::vector<double> dosage;
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
    double ss = std::accumulate(dosage.begin(), dosage.end(), 0.0)/dosage.size()/2.0;
    F1.push_back(1.0 - ss);    
    F2.push_back(ss);    

    // Eigen::VectorXd à partir du vecteur de dosages
    Eigen::Map<Eigen::VectorXd> SNP(&dosage[0], dosage.size());

    v = (PP.selfadjointView<Eigen::Lower>()*SNP).dot(SNP);
    t = SNP.dot(Py.col(0));
    s(i-beg) = t*t/v;

    dosage.clear();
  }
  
  List S;
  S["id"] = wrap(SNP_ID);
  S["chr"] = wrap(CHR);
  S["pos"] = wrap(POS);
  S["A1"] = wrap(AL1);
  S["A2"] = wrap(AL2);
  S["freq1"] = wrap(F1);
  S["freq2"] = wrap(F2);
  S["score"] = s;

  return S;
}


// float version
//[[Rcpp::export]]
List GWAS_dosage_lmm_score_f(CharacterVector filename, NumericVector PY, NumericMatrix P, int beg, int end) {

  dosages in(filename);

  int n = PY.size();
  if(P.nrow() != n || P.ncol() != n) 
    stop("Dimensions mismatch\n");

  Eigen::VectorXf Py(n);
  for(int i = 0; i < n; i++) Py(i) = PY[i];

  Eigen::MatrixXf PP(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) 
      PP(i,j) = P(i,j);
     

  int r = end-beg+1;

  Eigen::VectorXf SNP(n);
  NumericVector s(r);
  double t, v;  
  
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
    double ss = std::accumulate(dosage.begin(), dosage.end(), 0.0)/dosage.size()/2.0;
    F1.push_back(1.0 - ss);    
    F2.push_back(ss);    

    // Eigen::VectorXd à partir du vecteur de dosages
    Eigen::Map<Eigen::VectorXf> SNP(&dosage[0], dosage.size());

    v = (PP.selfadjointView<Eigen::Lower>()*SNP).dot(SNP);
    t = SNP.dot(Py);
    s(i-beg) = t*t/v;    

    dosage.clear();
  }
  
  List S;
  S["id"] = wrap(SNP_ID);
  S["chr"] = wrap(CHR);
  S["pos"] = wrap(POS);
  S["A1"] = wrap(AL1);
  S["A2"] = wrap(AL2);
  S["freq1"] = wrap(F1);
  S["freq2"] = wrap(F2);
  S["score"] = s;

  return S;
}


RcppExport SEXP GWAS_dosage_lmm_score(SEXP filenameSEXP, SEXP PYSEXP, SEXP PSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_dosage_lmm_score(filename, PY, P, beg, end));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP GWAS_dosage_lmm_score_f(SEXP filenameSEXP, SEXP PYSEXP, SEXP PSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PY(PYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_dosage_lmm_score_f(filename, PY, P, beg, end));
    return rcpp_result_gen;
END_RCPP
}

