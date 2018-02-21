#include "snp_filter.h"

using namespace Rcpp;

bool snp_filter::operator()(int chr, std::string & snp, int bp, double cm) {
  if(!H.m) return true; // empty hash
  if(H.htype == chr_pos) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
  if(H.htype == snpid) {
    int a = H.lookup(snp);
    return (a != NA_INTEGER);
  }
  return false; // ???
}

bool snp_filter::operator()(int chr, int bp) {
  int a = H.lookup(chr, bp);
  return (a != NA_INTEGER);
}

bool snp_filter::operator()(std::string & snp) {
  int a = H.lookup(snp);
  return (a != NA_INTEGER);
}

