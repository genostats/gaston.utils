#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "dosage_files.h"

using namespace Rcpp;

//[[Rcpp::export]]
void extract_dosage_file(std::string filename, std::string outname, LogicalVector keep) {
  dosages in( filename );
  ogzstream out( outname.c_str() );

  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> CHR, SNP_ID, AL1, AL2;
  std::vector<int> POS;
  std::vector<double> dosage;
  int nb_ind = keep.size();
  char buffer[10];


  // header line
  out << "id chr pos A1 A2 ";
  for(int i = 0; i < nb_ind; i++) {
    if(keep[i])
      out << " " << in.samples[i]; 
  }
  out << "\n";

  while( in.read_line(dosage, snp_id, snp_pos, chr, A1, A2) ) {
    // check for right number of dosages
    int dos_le = dosage.size();
    if(nb_ind != dos_le) {
      Rcerr << "While reading SNP #" << POS.size()+1 << " with id = " << snp_id << "\n";
      Rcerr << "Read " << dos_le << " dosages, instead of " << nb_ind << "\n";
      stop("Dimensions mismatch");
    }
    out << snp_id << " " << chr << " " << snp_pos << " " << A1 << " " << A2;
    for(int i = 0; i < nb_ind; i++) {
      if(keep[i]) {
        sprintf(buffer, "%.4f", dosage[i]);
        out << " " << buffer; 
      }
    }
    out << "\n";
    dosage.clear();
  }
  out.close();
}


RcppExport SEXP extract_dosage_file(SEXP filenameSEXP, SEXP outnameSEXP, SEXP keepSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type outname(outnameSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type keep(keepSEXP);
    extract_dosage_file(filename, outname, keep);
    return R_NilValue;
END_RCPP
}

