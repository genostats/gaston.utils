#include <Rcpp.h>

#ifndef CHR_CONVERT
#define CHR_CONVERT
using namespace Rcpp;

// UTILISER tolower() ?

static List chr_ids;

void set_chr_ids(List L);

inline int chr_to_int(std::string & chr) {
  int r = std::atoi(chr.c_str());
  if( r == 0 && chr_ids.containsElementNamed(chr.c_str()) )
    r = chr_ids[chr];
  return r;
}

inline int chr_to_int(char * chr) {
  int r = std::atoi(chr);
  if( r == 0 && chr_ids.containsElementNamed(chr) )
    r = chr_ids[chr];
  return r;
}



#endif
