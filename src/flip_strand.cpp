#include <Rcpp.h>
#include <string>
// #include "gaston/flip_strand.h"
#include "flip_strand.h"

using namespace Rcpp;


CharacterVector flip_strand(CharacterVector Allele) {
  std::vector<std::string> R;
  for(auto x : Allele) {
    R.push_back( flip_strand( CHAR(x) );
  }
  return wrap(R);
}
