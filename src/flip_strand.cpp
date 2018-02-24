#include <Rcpp.h>
#include <string>
// #include "gaston/flip_strand.h"
#include "flip_strand.h"

using namespace Rcpp;

//[[Rcpp::export]]
CharacterVector flip_strand(CharacterVector Allele) {
  std::vector<std::string> R;
  for(auto x : Allele) {
    R.push_back( flip_strand( CHAR(x) ) );
  }
  return wrap(R);
}

RcppExport SEXP gg_flip_strand(SEXP AlleleSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< CharacterVector >::type Allele(AlleleSEXP);
    return flip_strand(Allele);
END_RCPP
}

