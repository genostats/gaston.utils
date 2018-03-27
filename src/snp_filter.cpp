#include "snp_filter.h"
#include <string>

using namespace Rcpp;

// ------------------------------------------------------------------
bool snp_filter::operator()(std::string & snp, int chr, int bp) {

  if(t == nofilter || t == range_cm) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp)
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
 

  // if(!H.m) return true; // empty hash

  if(H.htype == chr_pos) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
  if(H.htype == snpid) {
    int a = H.lookup(snp);
    return (a != NA_INTEGER);
  }

  if(H.htype == snpid_chr_pos) {
    int a = H.lookup(snp, chr, bp);
    return (a != NA_INTEGER);
  }

  stop("Wrong hash type !");
  return false; // ???
}

// ------------------------------------------------------------------
bool snp_filter::operator()(std::string & snp, int chr, int bp, double cm) {

  if(t == nofilter) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp)
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
 
  if(t == range_cm) 
    return( chr == chr_ && low_cm <= cm && cm <= high_cm );
 

  // if(!H.m) return true; // empty hash

  if(H.htype == chr_pos) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
  if(H.htype == snpid) {
    int a = H.lookup(snp);
    return (a != NA_INTEGER);
  }

  if(H.htype == snpid_chr_pos) {
    int a = H.lookup(snp, chr, bp);
    return (a != NA_INTEGER);
  }

  stop("Wrong hash type !");
  return false; // ???
}

// ------------------------------------------------------------------
bool snp_filter::operator()(int chr, int bp) {
  if(t == nofilter) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp) 
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
  
  
  if(t != hash) return true;

  int a = H.lookup(chr, bp);
  return (a != NA_INTEGER);
}

// ------------------------------------------------------------------
bool snp_filter::operator()(std::string & snp) {
 
  if(t != hash) return true;

  int a = H.lookup(snp);
  return (a != NA_INTEGER);
}

