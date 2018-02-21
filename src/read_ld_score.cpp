#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "gzstream.h"
#include "snp_filter.h"

using namespace Rcpp;

List read_ld_score_filtered(std::vector<std::string> Filenames, snp_filter & F) {
  std::vector<int> CHR, BP;
  std::vector<std::string> SNP;
  std::vector<double> CM, MAF, L2;
 
  int total_nb_snps;

  for(std::string file : Filenames) {
    // open file
    igzstream in( file.c_str() );
    if(!in.good()) stop("Can't open file "+file);
  
    std::string line, str;
    // read and check header
    if(!std::getline(in, line))
      stop("File "+file+" is empty");
    
    std::vector<std::string> splitted;
    std::istringstream li(line);
    while(li >> str) splitted.push_back(str);
    if(splitted.size() != 6 || splitted[0] != "CHR" 
                            || splitted[1] != "SNP" 
                            || splitted[2] != "BP" 
                            || splitted[3] != "CM"
                            || splitted[4] != "MAF"
                            || splitted[5] != "L2")
       stop("Header of file "+file+" isn't what I expected");
  
     // let's go
     int chr, bp;
     std::string snp;
     double cm, maf, l2;
  
     while(std::getline(in, line)) {
       total_nb_snps++;
       li = std::istringstream(line);
       if(!(li >> chr >> snp >> bp >> cm >> maf >> l2)) {
         Rcerr << "last read line : [" << line << "]\n"; 
         stop("Format error in file "+file);
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
   }
   List D;
   D["CHR"] = wrap(CHR);
   D["SNP"] = wrap(SNP);
   D["BP"]  = wrap(BP);
   D["CM"]  = wrap(CM);
   D["MAF"] = wrap(MAF);
   D["L2"]  = wrap(L2);
   D["total_nb_snps"] = total_nb_snps;

   return D;
}

RcppExport SEXP read_ld_score_filt_chr_pos(SEXP FilenamesSEXP, SEXP CHRSEXP, SEXP POSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type Filenames(FilenamesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CHR(CHRSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type POS(POSSEXP);
    snp_filter F(CHR, POS);
    rcpp_result_gen = Rcpp::wrap(read_ld_score_filtered(Filenames, F));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP read_ld_score_filt_id(SEXP FilenamesSEXP, SEXP IDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type Filenames(FilenamesSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type ID(IDSEXP);
    snp_filter F(ID);
    rcpp_result_gen = Rcpp::wrap(read_ld_score_filtered(Filenames, F));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP read_ld_score_all(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type filename(filenameSEXP);
    snp_filter all;
    rcpp_result_gen = Rcpp::wrap(read_ld_score_filtered(filename, all));
    return rcpp_result_gen;
END_RCPP
}

