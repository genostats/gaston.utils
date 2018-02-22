#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "gzstream.h"
#include "snp_filter.h"

using namespace Rcpp;

List read_snp_list_filtered(std::vector<std::string> Filenames, snp_filter & F) {
  std::vector<std::string> SNP, A1, A2;
  
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
    if(splitted.size() != 3 || splitted[0] != "SNP" 
                            || splitted[1] != "A1" 
                            || splitted[2] != "A2")
       stop("Header of file "+file+" isn't what I expected");
  
     // let's go
     std::string snp, a1, a2;
  
     while(std::getline(in, line)) {
       std::istringstream li(line);
       if(!(li >> snp >> a1 >> a2)) {
         Rcerr << "last read line : [" << line << "]\n"; 
         stop("Format error in file "+file);
       }
       if(F(snp)) {
         SNP.push_back(snp);
         A1.push_back(a1);
         A2.push_back(a2);
       }
     }
   }
   List D;
   D["SNP"] = wrap(SNP);
   D["A1"]  = wrap(A1);
   D["A2"]  = wrap(A2);

   return D;
}

RcppExport SEXP read_snp_list_filt_id(SEXP FilenamesSEXP, SEXP IDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type Filenames(FilenamesSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type ID(IDSEXP);
    snp_filter F(ID);
    rcpp_result_gen = Rcpp::wrap(read_snp_list_filtered(Filenames, F));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP read_snp_list_all(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type filename(filenameSEXP);
    snp_filter all;
    rcpp_result_gen = Rcpp::wrap(read_snp_list_filtered(filename, all));
    return rcpp_result_gen;
END_RCPP
}

