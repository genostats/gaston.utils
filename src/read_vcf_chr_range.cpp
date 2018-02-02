#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include "matrix4.h"
#include <fstream>
#include "gzstream.h"
#include "read_vcf_line.h"
#include "token.h"
#include "default_value.h"

using namespace Rcpp;

void set(uint8_t * data, size_t j, uint8_t val) {
  uint8_t & a = data[j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  a |= (val << ((j%4)*2)); // set to val
}

void read_vcf_header(igzstream & in, std::vector<std::string> & samples, std::vector<std::string> & format_ids, std::vector<std::string> & info_ids) {
  // skip header and read samples
  std::string line, id;
  if(!std::getline(in, line))
    stop("File is empty");

  while(std::getline(in, line)) {
    if(line.substr(0,1) != "#") stop("Bad VCF format");
    if(line.substr(0,2) == "##") {
      if(line.substr(0,11) == "##INFO=<ID=") {
        std::istringstream li(line);
        std::getline(li, id, ',');
        info_ids.push_back(id.substr(11));
      }
      if(line.substr(0,13) == "##FORMAT=<ID=") {
        std::istringstream li(line);
        std::getline(li, id, ',');
        format_ids.push_back(id.substr(13));
      }
    } else {
      read_vcf_samples(line, samples);
      break; // fin
    }
  }
}

template<typename T>
void do_get_info(std::string & info, std::string & key, T & value) {
  Rcout << "info = " << info << "\n";
  Rcout << "key = " << key << "\n";
  size_t start = info.find(key+"=");
  if(start == info.npos) return;
  std::string s = info.substr(start + key.length() + 1);
  Rcout << "s = " << s << "\n";
  std::istringstream is(s);
  std::string val;
  std::getline(is, val, ';');
  Rcout << "val = " << val << "\n";
  value = sto<T>(val);
}

void parse_vcf_line_genotypes(std::string line, std::vector<std::string> & genotypes, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2, std::string & qual,
                     std::string & filter, std::string & info) {
  std::istringstream li(line);
  std::string format;
  if(!(li >> chr >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    Rcerr << line;
    stop("VCF file format error");
  }
  
  int pos = token_position(format, "GT");
  if(pos == 120) stop("bite");
  std::string G;
  while(li >> G) {
    genotypes.push_back( token_at_position<std::string>(G, pos) );
  }
}

// conversion d'un génotype 0|0 ou 0/0 -> 0 etc
inline int geno_conv(std::string s) {
  int g = 0;
  int le = s.length();
  if(le == 3) { // cas diploide
    if(s[0] == '1') g++;
    if(s[2] == '1') g++;
    if(s[0] == '.' || s[2] == '.') g = 3; // missing value : NA
  } else if(le == 1) { // cas haploide
    if(s[0] == '1') g++;
    if(s[0] == '.') g = 3; // missing value : NA
  } else { 
    g = 3;
  }
  return g;
} 

//[[Rcpp::export]]
List read_vcf_chr_range(std::string filename, bool get_info, std::string chr0, int low, int high) {
  // open file
  igzstream in( filename.c_str() );
  if(!in.good()) stop("Can't open file");

  // read header
  std::vector<std::string> samples_, format_ids, info_ids;
  read_vcf_header(in, samples_, format_ids, info_ids);

  // for(auto i: format_ids) Rcout << i << "\n";
  // for(auto i: info_ids) Rcout << i << "\n";

  CharacterVector samples = wrap(samples_);
  int nsamples = samples.length();

  List L;
  std::vector<std::string> chr, id, ref, alt, filter, qual;
  std::vector<int> pos;
  std::vector<std::vector<float>> INFO(info_ids.size());

  //  std::vector<double> qual;

  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé

  std::vector<uint8_t *> data;
  uint8_t * data_;

  int i = 0;
  std::string line;
  while(std::getline(in, line)) {
    int pos_ = 0;
    std::string chr_, id_("(no SNP read yet)"), ref_, alt_, filter_, info_, qual_;
    // double qual_;

    std::vector<std::string> genotypes;
    parse_vcf_line_genotypes(line, genotypes, id_, pos_, chr_, ref_, alt_, qual_, filter_, info_);
 
    if(chr0 != "") { // there's a chr condition
      Rcout << "(" << chr0 << ")(" << chr_ << ") -> " << chr0.compare(chr_) << "\n";
      if(chr0.compare(chr_)) // -> if this is different
        continue; // skip
      // note : high < 0 implies no range condition
      // take care of that in R code
      Rcout << "pos_ = " << pos_ << "\n";
      if(high > 0 && (pos_ < low || pos_ > high))
        continue; // skip
    }

    if(genotypes.size() != nsamples)
      Rf_error("VCF format error while reading SNP %s chr = %s pos %d", id_.c_str(), chr_.c_str(), pos_);

    // push back
    chr.push_back(chr_);
    pos.push_back(pos_);
    id.push_back(id_);
    ref.push_back(ref_);
    alt.push_back(alt_);
    qual.push_back(qual_); 
    filter.push_back(filter_);

    if(get_info) {
      for(int k = 0; k < info_ids.size(); k++) {
        float v = default_value<float>();
        do_get_info<float>(info_, info_ids[k], v);
        INFO[k].push_back(v);
      }
    }

    // nouvelle ligne de données    
    data_ = new uint8_t [pX->true_ncol];
    std::fill(data_, data_ + pX->true_ncol, 255); // c'est important de remplir avec 3 -> NA

    int j = 0;
    for(int j1 = 0; j1 < nsamples; j1++) {
      int g = geno_conv(genotypes[j1]);
      set(data_, j, g);
      j++;
    }
    data.push_back(data_);
    i++;
  }
  
  // et on finit la construction de la matrice
  pX->nrow = i; 
  if(i > 0) {
    pX->data = new uint8_t * [i];
    for(int j = 0; j < i; j++) pX->data[j] = data[j];
  }

  L["id"] = id;
  L["pos"] = pos;
  L["chr"] = chr;
  L["A1"] = ref;
  L["A2"] = alt;
  L["quality"] = qual;
  L["filter"] = filter;
  if(get_info) {
    for(int k = 0; k < info_ids.size(); k++) {
      L[ "info."+info_ids[k] ] = INFO[k];
    }
  }
  L["bed"] = pX;
  L["samples"] = samples;
  return L;
}

RcppExport SEXP gg_read_vcf_chr_range(SEXP filenameSEXP, SEXP get_infoSEXP, SEXP chr0SEXP, SEXP lowSEXP, SEXP highSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type get_info(get_infoSEXP);
    Rcpp::traits::input_parameter< std::string >::type chr0(chr0SEXP);
    Rcpp::traits::input_parameter< int >::type low(lowSEXP);
    Rcpp::traits::input_parameter< int >::type high(highSEXP);
    rcpp_result_gen = Rcpp::wrap(read_vcf_chr_range(filename, get_info, chr0, low, high));
    return rcpp_result_gen;
END_RCPP
}
