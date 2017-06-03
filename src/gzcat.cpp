#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"

using namespace Rcpp;

//[[Rcpp::export]]
void gzcat(CharacterVector filename) {
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");
  for(std::string line; std::getline(in, line); ) {
    Rcout << line << std::endl;
  }
  return;
}

RcppExport SEXP zz_gzcat(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    gzcat(filename);
    return R_NilValue;
END_RCPP
}


