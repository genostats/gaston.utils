#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "read_gen_file.h"

using namespace Rcpp;

List read_gen_file(CharacterVector filename) {
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");
  std::string snp_id, A1, A2;
  int snp_pos;
  std::vector<std::string> SNP_ID, AL1, AL2;
  std::vector<int> POS;
  std::vector<double> dosage;
  while( read_gen_line(in, dosage, snp_id, snp_pos, A1, A2) ) {
    SNP_ID.push_back(snp_id);
    POS.push_back(snp_pos);
    AL1.push_back(A1);
    AL2.push_back(A2);   
  }
  in.close();
  List L;
  L["snp.id"] = wrap(SNP_ID);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  NumericVector dos = wrap(dosage);
  dos.attr("dim") = Dimension( dosage.size()/POS.size(), POS.size() );
  L["dosages"] = dos;
  return L;
}

NumericVector gen_file_dim(CharacterVector filename) {
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");
  std::string snp_id, A1, A2;
  int snp_pos;
  std::vector<double> dosage;
  read_gen_line(in, dosage, snp_id, snp_pos, A1, A2);
  int nb_inds = dosage.size();
  int nb_snps = 1;
  dosage.clear();
  while( read_gen_line(in, dosage, snp_id, snp_pos, A1, A2) ) {
    nb_snps++;
    if(nb_inds != dosage.size())
      stop("gen file format error");
    dosage.clear();
  }
  in.close();
  return NumericVector::create(nb_inds, nb_snps);
}

//[[Rcpp::export]]
int nb_inds_gen_file(CharacterVector filename) {
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");
  std::string snp_id, A1, A2;
  int snp_pos;
  std::vector<double> dosage;
  read_gen_line(in, dosage, snp_id, snp_pos, A1, A2);
  return dosage.size();
}

RcppExport SEXP zz_read_gen_file(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_gen_file(filename));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP zz_gen_file_dim(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_file_dim(filename));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP zz_nb_inds_gen_file(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(nb_inds_gen_file(filename));
    return rcpp_result_gen;
END_RCPP
}

