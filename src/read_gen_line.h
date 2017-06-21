#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#ifndef GASTONread_gen_line
#define GASTONread_gen_line
using namespace Rcpp;

template<typename scalar_t>
inline void parse_gen_line(std::string line, std::vector<scalar_t> & dosage, std::string & snp_id, 
                          int & snp_pos, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  std::string tirets;
  if(!(li >> tirets && li >> snp_id && li >> snp_pos && li >> A1 && li >> A2)) 
    stop("gen file format error");
  scalar_t p0, p1, p2;
  while(li >> p0 && li >> p1 && li >> p2) {
    dosage.push_back(p1 + 2*p2);
  }
}

template<typename scalar_t>
inline bool read_gen_line(igzstream & in, std::vector<scalar_t> & dosage, std::string & snp_id, 
                          int & snp_pos, std::string & A1, std::string & A2) {
  std::string line;
  if(!std::getline(in, line)) {
    return false;
  }
  parse_gen_line(line, dosage, snp_id, snp_pos, A1, A2);
  return true;
}
#endif
