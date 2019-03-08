#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "token.h"
#include "stringstream_lite.h"
#include "chr_convert.h"
#include "snp_filter.h"
#include "flip_strand.h"

#ifndef GASTONread_vcf_line_allele_depth
#define GASTONread_vcf_line_allele_depth
using namespace Rcpp;

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

bool parse_vcf_line_allele_depth_filtered(std::string line, std::vector<int> & depth1, std::vector<int> & depth2, 
                     std::string & snp_id, int & snp_pos, int & chr, 
                     std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info, snp_filter & FILTER) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, "AD");
  if(pos < 0) stop("VCF error (No 'AD' found)");

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
    // maintenant le champ AD est dans la chaine pointée par tok.token, de longueur tok.token_length
    stringstream_lite tokAD(tok.token, ',');
    int D1, D2;
    if( !(tokAD >> D1 >> D2) )
      stop("VCF file format error (AD field ill formed)");

    if(swap) {
      depth1.push_back(D2);
      depth2.push_back(D1);
    } else {
      depth1.push_back(D1);
      depth2.push_back(D2);
    }
  }
  return true;
}

/*
// version 'filtered' + vecteur de booleens pour déterminer quels individus on prend
template<typename scalar_t>
bool parse_vcf_line_genotypes_filtered(std::string line, std::vector<scalar_t> & genotypes, 
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

  int pos = token_position(format, "AD");
  if(pos < 0) stop("VCF error (No 'AD' found)");

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
    // maintenant le champ AD est dans la chaine pointée par tok.token, de longueur tok.token_length
    if(which_samples[k++]) {
      stringstream_lite tokAD(tok.token, ',');
      int D1, D2;
      if( !(tokAD >> D1 >> D2) )
        stop("VCF file format error (AD field ill formed)");

      if(swap) {
        depth1.push_back(D2);
        depth2.push_back(D1);
      } else {
        depth1.push_back(D1);
        depth2.push_back(D2);
      }
    }
  }
  return true;
}
*/
#endif
