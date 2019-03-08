#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include "token.h"
#include "snp_filter.h"

#ifndef READ_FILTERED_
#define READ_FILTERED_

using namespace Rcpp;

void set(uint8_t * data, size_t j, uint8_t val);

template<typename T>
void do_get_info(std::string & info, std::string & key, T & value) {
  // Rcout << "info = " << info << "\n";
  // Rcout << "key = " << key << "\n";
  size_t start = info.find(key+"=");
  if(start == info.npos) return;
  std::string s = info.substr(start + key.length() + 1);
  // Rcout << "s = " << s << "\n";
  std::istringstream is(s);
  std::string val;
  std::getline(is, val, ';');
  // Rcout << "val = " << val << "\n";
  value = sto<T>(val);
}

bool all_equal(std::vector<std::string> a, std::vector<std::string> b);

void set_which_samples(std::vector<std::string> SAMPLES_, CharacterVector samples, 
       std::vector<bool> & which_samples, std::vector<std::string> & kept);

List read_vcf_filtered(std::vector<std::string> FILENAMES, bool get_info, snp_filter & FILTER, CharacterVector samples);

#endif
