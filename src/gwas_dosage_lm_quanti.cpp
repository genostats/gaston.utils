#include <Rcpp.h>
#include <RcppEigen.h>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"
#include <cmath>

using namespace Rcpp;

// Q = une matrice avec Q'Q = Id 
// --> décomposition QR de la matrice de covariables !
//[[Rcpp::export]]
List GWAS_dosage_lm_quanti(CharacterVector filename, NumericVector Y, NumericMatrix Q_, int beg, int end) {

  dosages in(filename);

  Eigen::Map<Eigen::MatrixXd> Q(as<Eigen::Map<Eigen::MatrixXd> >(Q_));
  Eigen::Map<Eigen::VectorXd> y(as<Eigen::Map<Eigen::VectorXd> >(Y));

  size_t n = Y.length();
  size_t p = Q.cols();
  // résidu Y par Q
  Eigen::VectorXd z(y - Q * (Q.transpose() * y));

  NumericVector beta(end-beg+1);
  NumericVector sd_beta(end-beg+1);

  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> SNP_ID, CHR, AL1, AL2;
  std::vector<int> POS;
  std::vector<double> dosage;

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

    // Eigen::VectorXd à partir du vecteur de dosages
    Eigen::Map<Eigen::VectorXd> SNP(&dosage[0], dosage.size());

    // résidu
    Eigen::VectorXd g(SNP - Q * (Q.transpose() * SNP));

    // régression de z sur g
    double gg = g.squaredNorm();
    double be = g.dot(z)/gg;
    beta(i-beg) = be;
    double s2 = (z - be*g).squaredNorm()/(n-p-1); // attention au p
    sd_beta(i-beg) = sqrt(s2/gg);
    dosage.clear();
  }
  List L;
  L["id"] = wrap(SNP_ID);
  L["chr"] = wrap(CHR);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  L["beta"] = beta;
  L["sd"] = sd_beta;
  return L;
}

RcppExport SEXP GWAS_dosage_lm_quanti(SEXP filenameSEXP, SEXP YSEXP, SEXP Q_SEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q_(Q_SEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_dosage_lm_quanti(filename, Y, Q_, beg, end));
    return rcpp_result_gen;
END_RCPP
}

