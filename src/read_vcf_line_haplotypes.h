/* 
 * 
 * idem read_vcf_filtered_haplos.cpp : bcp trop de copier coller là dedans 
 * puisque le seul truc qui change c'est geno_conv -> haplo_conv
 * 
 */

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "token.h"
#include "stringstream_lite.h"
#include "chr_convert.h"
#include "snp_filter.h"
#include "flip_strand.h"
#include "read_vcf_line.h"

#ifndef GASTONread_vcf_line_haplo
#define GASTONread_vcf_line_haplo
using namespace Rcpp;

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT


template<typename scalar_t>
inline std::pair<scalar_t, scalar_t> haplo_conv(char * s, int le) {
  scalar_t g1 = 0, g2 = 0;
  if(le == 3) { // cas diploide
    if(s[0] == '1') g1 = 1;
    if(s[2] == '1') g2 = 1;
    if(s[0] == '.') g1 = 3;
    if(s[2] == '.') g2 = 3; // missing value : NA
  } else if(le == 1) { // cas haploide
    g2 = 3; // pas de deuxième haplotype : NA
    if(s[0] == '1') g1 = 1;
    if(s[0] == '.') g1 = 3; // missing value : NA
  } else {
    g1 = g2 = 3;
  }
  return std::make_pair(g1, g2);
}





// version 'filtered'
template<typename scalar_t>
bool parse_vcf_line_haplotypes_filtered(std::string line, std::vector<scalar_t> & haplotypes, 
                     std::string & snp_id, int & snp_pos, int & chr, 
                     std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info, snp_filter & FILTER) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, "GT");
  if(pos < 0) stop("VCF error (No 'GT' found)");

  bool swap = false, flip = false;
  if(!FILTER(snp_id, chr, snp_pos, A1, A2, flip, swap)) 
    return false;

  if(swap) {
    std::string tmp(A1);
    A1 = A2;
    A2 = tmp;
  }
  if(flip) {
    A1 = flip_strand(A1);
    A2 = flip_strand(A2);
  }

  while( li.next_token() > 0 ) { // li.token pointe sur une chaîne avec le génotype en position pos
    stringstream_lite tok(li.token, ':'); // les champs sont séparés par des ':'
    for(int i = 0; i <= pos; i++) { // <= pos car même si pos = 0 il faut lire un token... 
      if(tok.next_token() == 0) 
        stop("VCF file format error");
    }

    // conversion du token t1 en génotype // ** tout ce bloc change
    auto g = haplo_conv<scalar_t>(tok.token, tok.token_length);
    if(swap) {
      haplotypes.push_back(1-g.first);
      haplotypes.push_back(1-g.second);
    } else {
      haplotypes.push_back(g.first);
      haplotypes.push_back(g.second);
    }
  }
  return true;
}


// version 'filtered' + vecteur de booleens pour déterminer quels individus on prend
template<typename scalar_t>
bool parse_vcf_line_haplotypes_filtered(std::string line, std::vector<scalar_t> & haplotypes, 
                     std::string & snp_id, int & snp_pos, int & chr, 
                     std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info, snp_filter & FILTER, 
                     std::vector<bool> & which_samples) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, "GT");
  if(pos < 0) stop("VCF error (No 'GT' found)");

  bool swap = false, flip = false;
  if(!FILTER(snp_id, chr, snp_pos, A1, A2, flip, swap)) 
    return false;

  if(swap) {
    std::string tmp(A1);
    A1 = A2;
    A2 = tmp;
  }
  if(flip) {
    A1 = flip_strand(A1);
    A2 = flip_strand(A2);
  }

  int k = 0;
  while( li.next_token() > 0 ) { // li.token pointe sur une chaîne avec le génotype en position pos
    stringstream_lite tok(li.token, ':'); // les champs sont séparés par des ':'
    for(int i = 0; i <= pos; i++) { // <= pos car même si pos = 0 il faut lire un token... 
      if(tok.next_token() == 0) 
        stop("VCF file format error");
    }

    if(which_samples[k++]) {
      // conversion du token t1 en génotype // ** tout ce bloc change
      auto g = haplo_conv<scalar_t>(tok.token, tok.token_length);
      if(swap) {
        haplotypes.push_back(1-g.first);
        haplotypes.push_back(1-g.second);
      } else {
        haplotypes.push_back(g.first);
        haplotypes.push_back(g.second);
      }
    }
  }
  return true;
}

#endif
