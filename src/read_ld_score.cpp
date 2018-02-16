#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "gzstream.h"
#include "gaston/snp_hash.h"

using namespace Rcpp;

class snp_filter {
  public:
  virtual bool operator()(int chr, std::string snp, int bp, double cm) {
    return true;
  }
};

class snp_filter_by_chr_pos : public snp_filter {
  SNPhash H;
  public:
  snp_filter_by_chr_pos(IntegerVector CHR, IntegerVector POS) : H(CHR, POS) { }
  bool operator()(int chr, std::string snp, int bp, double cm) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
};

List read_ld_score_filtered(std::string filename, snp_filter & F) {
  // open file
  igzstream in( filename.c_str() );
  if(!in.good()) stop("Can't open file");

  std::string line, str;
  // read and check header
  if(!std::getline(in, line))
    stop("File is empty");
  
  std::vector<std::string> splitted;
  std::istringstream li(line);
  while(li >> str) splitted.push_back(str);
  if(splitted.size() != 6 || splitted[0] != "CHR" 
                          || splitted[1] != "SNP" 
                          || splitted[2] != "BP" 
                          || splitted[3] != "CM"
                          || splitted[4] != "MAF"
                          || splitted[5] != "L2")
     stop("Header isn't what I expected\n");

   // let's go
   int chr, bp;
   std::string snp;
   double cm, maf, l2;

   std::vector<int> CHR, BP;
   std::vector<std::string> SNP;
   std::vector<double> CM, MAF, L2;

   while(std::getline(in, line)) {
     li = std::istringstream(line);
     if(!(li >> chr >> snp >> bp >> cm >> maf >> l2)) {
       Rcerr << "line : [" << line << "]\n"; 
       stop("Format error\n");
     }
     if(F(chr, snp, bp, cm)) {
       CHR.push_back(chr);
       SNP.push_back(snp);
       BP.push_back(bp);
       CM.push_back(cm);
       MAF.push_back(maf);
       L2.push_back(l2);   
     }
   }

   List D;
   D["CHR"] = wrap(CHR);
   D["SNP"] = wrap(SNP);
   D["BP"]  = wrap(BP);
   D["CM"]  = wrap(CM);
   D["MAF"] = wrap(MAF);
   D["L2"]  = wrap(L2);

   return D;
}

//[[Rcpp::export]]
List read_ld_score_filt(std::string filename, IntegerVector CHR, IntegerVector POS) {
  snp_filter_by_chr_pos F(CHR, POS);
  return read_ld_score_filtered(filename, F);
}

RcppExport SEXP read_ld_score_all(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    snp_filter all;
    rcpp_result_gen = Rcpp::wrap(read_ld_score_filtered(filename, all));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP read_ld_score_filt(SEXP filenameSEXP, SEXP CHRSEXP, SEXP POSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CHR(CHRSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type POS(POSSEXP);
    rcpp_result_gen = Rcpp::wrap(read_ld_score_filt(filename, CHR, POS));
    return rcpp_result_gen;
END_RCPP
}
