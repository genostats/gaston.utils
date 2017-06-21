#include <Rcpp.h>
#include <RcppEigen.h>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"
#include <cmath>
#include "diago2_full.h"
#include "diago2_full_nocovar.h"

using namespace Rcpp;

// laisser en double ça va aussi vite (plus vite ?) et ça fait vraiment
// une différence si il y a des covariables
#define scalar double
#define MATRIX MatrixXd
#define VECTOR VectorXd

//[[Rcpp::export]]
List GWAS_dosage_lmm_lrt(CharacterVector filename, NumericVector Y, NumericMatrix X, 
                  int p, NumericVector Sigma, NumericMatrix U, int beg, int end, double tol) {

  dosages in(filename);

  int n = Sigma.size();
  int r = X.ncol();
  int max_iter = 10;

  if(Y.size() != n || X.nrow() != n || U.nrow() != n || U.ncol() != n)
    stop("Dimensions mismatch");

  // conversion en float...
  MATRIX y0(n, 1);
  for(int i = 0; i < n; i++)
    y0(i,0) = Y[i];

  MATRIX x0(n, r);
  for(int j = 0; j < r; j++)
    for(int i = 0; i < n; i++)
      x0(i,j) = X(i,j);

  VECTOR sigma(n);
  for(int i = 0; i < n; i++)
    sigma[i] = Sigma[i];

  MATRIX u(n,n);
  for(int j = 0; j < n; j++)
    for(int i = 0; i < n; i++)
      u(i,j) = U(i,j);

  MATRIX x = u.transpose() * x0;
  MATRIX y = u.transpose() * y0;

  // Zecteur SNPs
  VECTOR SNP(n);

  // declare vectors containing result
  NumericVector H2(end-beg+1);
  NumericVector LRT(end-beg+1);

  scalar h2 = 0;
  scalar likelihood0;

  // on commence par le modèle null
  if(r == 1) { 
    // object for likelihood maximization
    diag_full_likelihood_nocovar<MATRIX, VECTOR, scalar> A(p, y, sigma);
    A.newton_max(h2, 0, 0.99, tol, max_iter, false);
    likelihood0 = A.likelihood();
  } else { 
    MATRIX x1 = x.leftCols(r-1);
    // object for likelihood maximization
    diag_full_likelihood<MATRIX, VECTOR, scalar> A(p, y, x1, sigma);
    A.newton_max(h2, 0, 0.99, tol, max_iter, false);
    likelihood0 = A.likelihood();
  }
  // object for likelihood maximization
  diag_full_likelihood<MATRIX, VECTOR, scalar> A(p, y, x, sigma);


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

    A.X.col(r-1) = u.transpose() * SNP;

    // likelihood maximization
    h2 = (h2 > 0.9)?0.9:h2;
    A.newton_max( h2, 0, 0.99, tol, max_iter, false);

    H2(i-beg) = h2;
    LRT(i-beg) = 2*(A.likelihood() - likelihood0);

    dosage.clear();
  }

  List R;
  R["id"] = wrap(SNP_ID);
  R["chr"] = wrap(CHR);
  R["pos"] = wrap(POS);
  R["A1"] = wrap(AL1);
  R["A2"] = wrap(AL2);
  R["h2"] = H2;
  R["LRT"] = LRT;
  return R;
}

RcppExport SEXP GWAS_dosage_lmm_lrt(SEXP filenameSEXP, SEXP YSEXP, SEXP XSEXP, SEXP pSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP begSEXP, SEXP endSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS_dosage_lmm_lrt(filename, Y, X, p, Sigma, U, beg, end, tol));
    return rcpp_result_gen;
END_RCPP
}

