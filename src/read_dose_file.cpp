#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"

using namespace Rcpp;

List read_dose_file(CharacterVector filename) {
  dosages in( filename );
  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> CHR, SNP_ID, AL1, AL2;
  std::vector<int> POS;
  std::vector<double> dosage;
  while( in.read_line(dosage, snp_id, snp_pos, chr, A1, A2) ) {
    SNP_ID.push_back(snp_id);
    POS.push_back(snp_pos);
    CHR.push_back(chr);
    AL1.push_back(A1);
    AL2.push_back(A2);   
  }
  List L;
  L["snp.id"] = wrap(SNP_ID);
  L["chr"] = wrap(CHR);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  if(in.type == 10) L["samples"] = wrap(in.samples); // VCF
  NumericVector dos = wrap(dosage);
  dos.attr("dim") = Dimension( dosage.size()/POS.size(), POS.size() );
  L["dosages"] = dos;
  return L;
}

NumericVector dose_file_dim(CharacterVector filename) {
  dosages in( filename );
  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<double> dosage;
  in.read_line(dosage, snp_id, snp_pos, chr, A1, A2);
  int nb_inds = dosage.size();
  int nb_snps = 1;
  dosage.clear();
  while( in.read_line(dosage, snp_id, snp_pos, chr, A1, A2) ) {
    nb_snps++;
    if(nb_inds != dosage.size())
      stop("File format error");
    dosage.clear();
  }
  return NumericVector::create(nb_inds, nb_snps);
}

int nb_inds_dose_file(CharacterVector filename) {
  dosages in( filename );
  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<double> dosage;
  in.read_line(dosage, snp_id, snp_pos, chr, A1, A2);
  return dosage.size();
}

RcppExport SEXP zz_read_dose_file(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_dose_file(filename));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP zz_dose_file_dim(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(dose_file_dim(filename));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP zz_nb_inds_dose_file(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(nb_inds_dose_file(filename));
    return rcpp_result_gen;
END_RCPP
}

