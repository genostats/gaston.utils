#include <Rcpp.h>
// #include "gaston/snp_hash.h"
#include "snp_hash.h"

using namespace Rcpp;

#ifndef SNP_FILTER
#define SNP_FILTER

class snp_filter {
  SNPhash H;
  public:

  snp_filter() {}
  snp_filter(IntegerVector CHR, IntegerVector POS) : H(CHR, POS) { }
  snp_filter(CharacterVector ID) : H(ID) { }
  bool operator()(int chr, std::string & snp, int bp, double cm);

  bool operator()(int chr, int bp);
  bool operator()(std::string & snp);
};

#endif
