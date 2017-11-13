#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "token.h"
#include "read_vcf_line.h"

void read_vcf_samples(std::string line, std::vector<std::string> & samples) {
  std::istringstream li(line);
  std::string G;
  for(int i = 0; i < 9; i++) {
    if(!(li >> G))
      stop("VCF file format error");
  }

  while(li >> G) {
    samples.push_back( G );
  }
}
