#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include "matrix4.h"
#include <fstream>
#include "gzstream.h"
#include "read_vcf_line.h"

using namespace Rcpp;

void set(uint8_t * data, size_t j, uint8_t val) {
  uint8_t & a = data[j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  a |= (val << ((j%4)*2)); // set to val
}

void parse_vcf_line_genotypes(std::string line, std::vector<std::string> & genotypes, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info) {
  std::istringstream li(line);
  std::string format;
  if(!(li >> chr >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format))
    stop("VCF file format error");
  
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
List read_vcf_chr_range(CharacterVector filename, bool get_info) {
  // open file
  if(filename.length() != 1) stop("filename should be a CharacterVector of length 1");
  igzstream in( (char *) filename[0] );
  if(!in.good()) stop("Can't open file");

  // skip header and read samples
  std::string line;
  if(!std::getline(in, line))
    stop("File is empty");

  std::vector<std::string> samples_;
  while(std::getline(in, line)) {
    if(line.substr(0,1) != "#") stop("Bad VCF format");
    if(line.substr(0,2) != "##") {
      read_vcf_samples(line, samples_);
      break; // fin
    }
  }

  CharacterVector samples = wrap(samples_);
  int nsamples = samples.length();

  List L;
  std::vector<std::string> chr, id, ref, alt, filter, info;
  std::vector<int> pos;
  std::vector<double> qual;

  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé

  std::vector<uint8_t *> data;
  uint8_t * data_;

  int i = 0;
  while(std::getline(in, line)) {
    int pos_ = 0;
    std::string chr_, id_("(no SNP read yet)"), ref_, alt_, filter_, info_;
    double qual_;

    std::vector<std::string> genotypes;
    parse_vcf_line_genotypes(line, genotypes, id_, pos_, chr_, ref_, alt_, qual_, filter_, info_);
      
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

    if(get_info) info.push_back(info_);
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
  if(get_info) L["info"] = info;

  L["bed"] = pX;
  L["samples"] = samples;
  return L;
}

