#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "ms_files.h"
#include "gaston/matrix4.h"

using namespace Rcpp;

//[[Rcpp::export]]
List read_ms_file(std::string filename) {
  msFile in(filename);
/*
  Rcout << "nsamples = " << in.nsamples << "\n";
  Rcout << "nreplicates = " << in.nreplicates << "\n";
  Rcout << "total segsites = " << in.total_segsites << "\n";

  Rcout << "rep = " << in.rep << "\n";
  Rcout << "segsites = " << in.segsites << "\n";
*/
  XPtr<matrix4> pX(new matrix4(in.total_segsites, in.nsamples));  

  std::vector<char> x;
  int i = 0;

  while( in.read_SNP(x) ) {
    // for(auto u : x) Rcout << (int) u;
    // Rcout << "\n";
    for(int j = 0; j < in.nsamples; j++)
      pX->set(i, j, x[j]);
    i++;
    if(i > in.total_segsites) 
      stop("ms file format error");
  }

  List L; 
  L["bed"] = pX;
  L["snps"] = in.total_segsites;
  L["reps"] = in.SEGSITES;
  L["haps"] = in.nsamples;

  return L;
}


RcppExport SEXP read_ms_file(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_ms_file(filename));
    return rcpp_result_gen;
END_RCPP
}

