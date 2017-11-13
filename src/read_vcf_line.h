#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "token.h"

#ifndef GASTONread_vcf_line
#define GASTONread_vcf_line
using namespace Rcpp;

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

template<typename scalar_t> 
void parse_vcf_line_dosages(std::string line, std::vector<scalar_t> & dosage, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  std::string qual, filter, info, format;
  if(!(li >> chr >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format))
    stop("VCF file format error");
  
  int pos = token_position(format, "DS");

  std::string G;
  while(li >> G) {
    dosage.push_back( token_at_position<scalar_t>(G, pos) );
  }
}


template<typename scalar_t> 
bool read_vcf_line_dosages(igzstream & in, std::vector<scalar_t> & dosage, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  std::string line;
  if(!std::getline(in, line))
    return false;
  parse_vcf_line_dosages<scalar_t>(line, dosage, snp_id, snp_pos, chr, A1, A2);
  return true;
}

void read_vcf_samples(std::string line, std::vector<std::string> & samples);
#endif
