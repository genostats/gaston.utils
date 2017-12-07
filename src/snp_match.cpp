// [[Rcpp::depends(BH)]]
//
#include <Rcpp.h>
#include <iostream>
#include "snp_hash.h"

using namespace Rcpp;

inline char flip_strand(char x) {
  if(x == 'A') return 'T';
  if(x == 'C') return 'G';
  if(x == 'G') return 'C';
  if(x == 'T') return 'A';
  return x;
}

std::string flip_strand(const char * str) {
  std::string s;
  while(*str) {
    s += flip_strand(*str++);
  }
  return s;
}

//[[Rcpp::export]]
IntegerVector which_duplicated_chr_pos(IntegerVector Chr1, IntegerVector Pos1) {
  SNPhash h(Chr1, Pos1);
  return wrap(h.dup_indices);
}

// snps_table est supposé avoir toutes les colonnes id chr pos a1 a2
// (en fait il suffit d'avoir les colonnes présentes dans x)
// On choisit le type de hash en fonction des colonnes présentes dans x
// (plus y en a plus on en met...)
// [[Rcpp::export]]
List SNPmatch(DataFrame x, DataFrame table) {
  hash_type ht;
  // determines hash type from x columns
  if(x.containsElementNamed("id")) {
    if(x.containsElementNamed("chr")) { 
      if(x.containsElementNamed("pos")) {
        if(x.containsElementNamed("A1") && x.containsElementNamed("A2")) {
          ht = snpid_chr_pos_al;
        } else {
           ht= snpid_chr_pos;
        }
      } else {
        ht = snpid;
      }
    } else { 
      ht = snpid;
    }
  } else {
    if(x.containsElementNamed("chr")) { 
      if(x.containsElementNamed("pos")) {
        if(x.containsElementNamed("A1") && x.containsElementNamed("A2")) {
          ht = chr_pos_al;
        } else {
          ht = chr_pos;
        }
      } else {
        stop("x can't be used for SNP matching");
      }
    } else { 
      stop("x can't be used for SNP matching");
    }
  }

  // creates hash from table columns
  // and run look up functions 
  List L;
  int n = x.nrows();
  SNPhash h;
  if(ht == snpid_chr_pos_al) { /****************************************/
    h = SNPhash( as<CharacterVector>(table["id"]), 
                 as<IntegerVector>(table["chr"]), 
                 as<IntegerVector>(table["pos"]), 
                 as<CharacterVector>(table["A1"]), 
                 as<CharacterVector>(table["A2"]));
    IntegerVector I(n);
    LogicalVector Swap(n);
    LogicalVector Flip(n);
    CharacterVector ID(as<CharacterVector>(x["id"]));
    IntegerVector CHR(as<IntegerVector>(x["chr"]));
    IntegerVector POS(as<IntegerVector>(x["pos"]));
    CharacterVector A1(as<CharacterVector>(x["A1"]));
    CharacterVector A2(as<CharacterVector>(x["A2"]));

    for(int i = 0; i < n; i++) {
      int ii = h.lookup( ID[i], CHR[i], POS[i], A1[i], A2[i] );
      if(ii != NA_INTEGER) {
        I[i] = ii;
        Swap[i] = false;
        Flip[i] = false;
        continue;
      }
      // try swapping alleles
      ii = h.lookup( ID[i], CHR[i], POS[i], A2[i], A1[i] );
      if(ii != NA_INTEGER) {
        I[i] = ii;
        Swap[i] = true;
        Flip[i] = false;
        continue;
      }
      // try flipping strand
      std::string b1 = flip_strand(CHAR(A1[i]));
      std::string b2 = flip_strand(CHAR(A2[i]));
      ii = h.lookup( ID[i], CHR[i], POS[i], b1, b2 );
      if(ii != NA_INTEGER) {
        I[i] = ii;
        Swap[i] = false;
        Flip[i] = true;
        continue;
      }
      // try flipping strand + swapping alleles
      ii = h.lookup( ID[i], CHR[i], POS[i], b2, b1);
      I[i] = ii;
      if(ii != NA_INTEGER) {
        Swap[i] = true;
        Flip[i] = true;
      }
    }
    // L["hash_index"] = wrap(h.index);
    L["index"] = I;
    L["swap"] = Swap;
    L["flip"] = Flip;
  } else if(ht == snpid_chr_pos) { /****************************************/
    h = SNPhash( as<CharacterVector>(table["id"]), 
                 as<IntegerVector>(table["chr"]), 
                 as<IntegerVector>(table["pos"]));
    IntegerVector I(n);
    CharacterVector ID(as<CharacterVector>(x["id"]));
    IntegerVector CHR(as<IntegerVector>(x["chr"]));
    IntegerVector POS(as<IntegerVector>(x["pos"]));
    for(int i = 0; i < n; i++) {
      I[i] = h.lookup( ID[i], CHR[i], POS[i]);
    }
    // L["hash_index"] = wrap(h.index);
    L["index"] = I;
  } else if(ht == snpid) {         /****************************************/
    h = SNPhash( as<CharacterVector>(table["id"])); 
    IntegerVector I(n);
    CharacterVector ID(as<CharacterVector>(x["id"]));
    for(int i = 0; i < n; i++) {
      I[i] = h.lookup( ID[i]);
    }
    // L["hash_index"] = wrap(h.index);
    L["index"] = I;
  } else if(ht == chr_pos_al) {   /****************************************/
    h = SNPhash( as<IntegerVector>(table["chr"]), 
                 as<IntegerVector>(table["pos"]), 
                 as<CharacterVector>(table["A1"]), 
                 as<CharacterVector>(table["A2"]));
    IntegerVector I(n);
    LogicalVector Swap(n);
    LogicalVector Flip(n);
    IntegerVector CHR(as<IntegerVector>(x["chr"]));
    IntegerVector POS(as<IntegerVector>(x["pos"]));
    CharacterVector A1(as<CharacterVector>(x["A1"]));
    CharacterVector A2(as<CharacterVector>(x["A2"]));

    for(int i = 0; i < n; i++) {
      int ii = h.lookup(CHR[i], POS[i], A1[i], A2[i] );
      if(ii != NA_INTEGER) {
        I[i] = ii;
        Swap[i] = false;
        Flip[i] = false;
        continue;
      }
      // try swapping alleles
      ii = h.lookup(CHR[i], POS[i], A2[i], A1[i] );
      if(ii != NA_INTEGER) {
        I[i] = ii;
        Swap[i] = true;
        Flip[i] = false;
        continue;
      }
      // try flipping strand
      std::string b1 = flip_strand(CHAR(A1[i]));
      std::string b2 = flip_strand(CHAR(A2[i]));
      ii = h.lookup(CHR[i], POS[i], b1, b2 );
      if(ii != NA_INTEGER) {
        I[i] = ii;
        Swap[i] = false;
        Flip[i] = true;
        continue;
      }
      // try flipping strand + swapping alleles
      ii = h.lookup(CHR[i], POS[i], b2, b1);
      I[i] = ii;
      if(ii != NA_INTEGER) {
        Swap[i] = true;
        Flip[i] = true;
      }
    }
    // L["hash_index"] = wrap(h.index);
    L["index"] = I;
    L["swap"] = Swap;
    L["flip"] = Flip;
  } else if(ht == chr_pos) {       /****************************************/
    h = SNPhash( as<IntegerVector>(table["chr"]), 
                 as<IntegerVector>(table["pos"]));
    IntegerVector I(n);
    IntegerVector CHR(as<IntegerVector>(x["chr"]));
    IntegerVector POS(as<IntegerVector>(x["pos"]));
    for(int i = 0; i < n; i++) {
      I[i] = h.lookup(CHR[i], POS[i]);
    }
    // L["hash_index"] = wrap(h.index);
    L["index"] = I;
  }
  return L;
}


RcppExport SEXP gg_which_duplicated_chr_pos(SEXP Chr1SEXP, SEXP Pos1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Chr1(Chr1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Pos1(Pos1SEXP);
    rcpp_result_gen = Rcpp::wrap(which_duplicated_chr_pos(Chr1, Pos1));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_SNPmatch(SEXP xSEXP, SEXP tableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type table(tableSEXP);
    rcpp_result_gen = Rcpp::wrap(SNPmatch(x, table));
    return rcpp_result_gen;
END_RCPP
}

