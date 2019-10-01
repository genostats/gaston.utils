#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "ms_files.h"
#include "gaston/matrix4.h"

using namespace Rcpp;

//[[Rcpp::export]]
List read_ms_file_haplotype(std::string filename, int rep) {
  msFile in(filename, rep);
  // Rcout << "nsamples = " << in.nsamples << "\n";
  // Rcout << "segsites = " << in.segsites << "\n";
  // for(auto p : in.pos) Rcout << p << " ";
  // Rcout << "\n";

  XPtr<matrix4> pX(new matrix4(in.segsites, in.nsamples));  

  std::vector<char> x;
  int j = 0;

  while( in.read_haplotype(x) ) {
    // for(auto u : x) Rcout << (int) u;
    // Rcout << "\n";
    for(int i = 0; i < in.segsites; i++)
      pX->set(i, j, x[i]);
    j++;
    if(j > in.nsamples) 
      stop("Too much haplotypes in ms file");
  }
  if(j < in.nsamples) 
    stop("Too few haplotypes in ms file");

  List L;
  L["bed"] = pX;
  L["snps"] = in.segsites;
  L["haps"] = in.nsamples;
  L["pos"] = in.pos;

  return L;
}


RcppExport SEXP read_ms_file_haplotype(SEXP filenameSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(read_ms_file_haplotype(filename, rep));
    return rcpp_result_gen;
END_RCPP
}

