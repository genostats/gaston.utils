#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include "gaston/matrix4.h"
#include <fstream>
#include "gzstream.h"
#include "read_vcf_line.h"
#include "token.h"
#include "default_value.h"
#include "snp_filter.h"
#include "read_vcf_header.h"
#include "read_vcf_filtered.h"

#define SHOW(a) Rcout << #a << " = " << a << "\n";

using namespace Rcpp;

// un set "sécurisé"
void set(uint8_t * data, size_t j, uint8_t val) {
  uint8_t & a = data[j/4];
  a &= ~(3 << ((j%4)*2));  // set to 00
  switch(val) {
    case 0:
      return;
    case 1:
      a |= (1 << ((j%4)*2)); // set to 1
      return;
    case 2:
      a |= (2 << ((j%4)*2)); // set to 2
      return;
    default:
      a |= (3 << ((j%4)*2)); // set to 3 = NA
  }
}

bool all_equal(std::vector<std::string> a, std::vector<std::string> b) {
  int n = a.size();
  if(n != b.size())
    return false;
  for(int i = 0; i < n; i++) 
    if( a[i].compare(b[i]) != 0 ) 
      return false;
  return true;
}

// on passe par des CharacterVector pour utiliser IndexHash...
// mais pour la suite je prends des std::vector... (que d'aller-retour)
// SAMPLES_ = les échantillons présents dans le fichier
// samples = les échantillons à garder
// OUTPUT : which_samples , booléens de logueur SAMPLES_ qui donne s'il faut ou non garder l'échantillon
//        : kept = le nom des échantillons gardés, dans l'ordre où ils apparaissent dans SAMPLES_
void set_which_samples(std::vector<std::string> SAMPLES_, CharacterVector samples, std::vector<bool> & which_samples, std::vector<std::string> & kept) {
  sugar::IndexHash<STRSXP> H(samples);
  H.fill();
  which_samples.resize(0);
  which_samples.reserve(SAMPLES_.size());

  kept.resize(0);
  kept.reserve(SAMPLES_.size());

  CharacterVector SAMPLES(wrap(SAMPLES_));
  for(auto x : SAMPLES) {
     if(H.contains(x)) {
       which_samples.push_back(true);
       kept.push_back(as<std::string>(x));
     } else
       which_samples.push_back(false);
  }

}

//[[Rcpp::export]]
List read_vcf_filtered(std::vector<std::string> FILENAMES, bool get_info, snp_filter & FILTER, CharacterVector samples) {

  if(FILENAMES.size() < 1)
    stop("Empty Filenames vector");

  std::vector<std::string> id, ref, alt, filter;
  std::vector<int> pos, chr;
  std::vector<double> qual;

  std::vector<std::string> SAMPLES, SAMPLES_kept, FORMAT_IDS, INFO_IDS;
  std::vector<bool> which_samples;
  bool filter_sample = (samples.size() > 0);

  // ***** begin to read first file ... ***************
  {
    std::string filename = FILENAMES[0];

    // open file
    igzstream in( filename.c_str() );
    if(!in.good()) stop("Can't open file "+filename);

    // read header
    read_vcf_header(in, SAMPLES, FORMAT_IDS, INFO_IDS);
  }
  // ***** we got what we need ! **********************

  // we can create the bed matrix and restart

  // échantillons à lire...
  if(filter_sample) {
    set_which_samples(SAMPLES, samples, which_samples, SAMPLES_kept);
  } else {
    SAMPLES_kept = SAMPLES;
  }

  int nsamples = SAMPLES_kept.size();

  std::vector<std::vector<float>> INFO( INFO_IDS.size() );

  XPtr<matrix4> pX(new matrix4(0, nsamples));  // avec nrow = 0 allocations() n'est pas appelé
  std::vector<uint8_t *> data;
  int nb_snps = 0;

  for(std::string filename : FILENAMES) {

    // open file
    igzstream in( filename.c_str() );
    if(!in.good()) stop("Can't open file "+filename);

    // read header
    std::vector<std::string> samples_, samples_kept, format_ids, info_ids;
    read_vcf_header(in, samples_, format_ids, info_ids);

    // échantillons à lire...
    if(filter_sample) {
      set_which_samples(samples_, samples, which_samples, samples_kept);
    } else {
      samples_kept = samples_;
    }

    // check if samples are the same ...
    if(!all_equal(SAMPLES_kept, samples_kept))
      stop("Samples should be present and in the same order in all VCF files");

    // et si on veut les info...
    if(get_info && !all_equal(INFO_IDS, info_ids))
      stop("With get_info = TRUE, the info IDs should be the same in all VCF files");

    uint8_t * data_;

    std::string line;
    while(std::getline(in, line)) {
      int pos_ = 0;
      std::string id_("(no SNP read yet)"), ref_, alt_, filter_, info_;
      int chr_;
      double qual_;

      std::vector<int> genotypes;
      if(filter_sample) {
        if(!parse_vcf_line_genotypes_filtered(line, genotypes, id_, pos_, chr_, ref_, alt_, qual_, filter_, info_, FILTER, which_samples))
          continue; // skip
      } else {
        if(!parse_vcf_line_genotypes_filtered(line, genotypes, id_, pos_, chr_, ref_, alt_, qual_, filter_, info_, FILTER))
          continue; // skip
      }

      if(genotypes.size() != nsamples)
        Rf_error("VCF format error while reading SNP %s chr = %d pos %d", id_.c_str(), chr_, pos_);

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
        int g = genotypes[j1];
        set(data_, j, g);
        j++;
      }
      data.push_back(data_);
      nb_snps++;
    }
    in.close();
  }
  
  // et on finit la construction de la matrice
  pX->nrow = nb_snps; 
  if(nb_snps > 0) {
    pX->data = new uint8_t * [nb_snps];
    for(int j = 0; j < nb_snps; j++) pX->data[j] = data[j];
  }

  List L;
  L["id"] = id;
  L["pos"] = pos;
  L["chr"] = chr;
  L["A1"] = ref;
  L["A2"] = alt;
  L["quality"] = qual;
  L["filter"] = filter;
  if(get_info) {
    for(int k = 0; k < INFO_IDS.size(); k++) {
      L[ "info."+INFO_IDS[k] ] = INFO[k];
    }
  }
  L["bed"] = pX;
  L["samples"] = SAMPLES_kept;
  return L;
}


RcppExport SEXP gg_read_vcf_chr_range(SEXP filenameSEXP, SEXP get_infoSEXP, SEXP chr0SEXP, SEXP lowSEXP, SEXP highSEXP, SEXP samples_SXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type get_info(get_infoSEXP);
    Rcpp::traits::input_parameter< int >::type chr(chr0SEXP);
    Rcpp::traits::input_parameter< int >::type low(lowSEXP);
    Rcpp::traits::input_parameter< int >::type high(highSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type samples(samples_SXP);
    if(chr < 0) {
      snp_filter F;
      rcpp_result_gen = Rcpp::wrap(read_vcf_filtered(filename, get_info, F, samples));
    } else if(high < 0) {
      snp_filter F(chr);
      rcpp_result_gen = Rcpp::wrap(read_vcf_filtered(filename, get_info, F, samples));
    } else {
      snp_filter F(chr, low, high);
      rcpp_result_gen = Rcpp::wrap(read_vcf_filtered(filename, get_info, F, samples));
    }
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_read_vcf_chr_pos(SEXP filenameSEXP, SEXP get_infoSEXP, SEXP chrSXP, SEXP posSXP, SEXP samples_SXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type get_info(get_infoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pos(posSXP);
    Rcpp::traits::input_parameter< CharacterVector >::type samples(samples_SXP);
    snp_filter F(chr, pos);
    rcpp_result_gen = Rcpp::wrap(read_vcf_filtered(filename, get_info, F, samples));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP gg_read_vcf_chr_pos_al(SEXP filenameSEXP, SEXP get_infoSEXP, SEXP chrSXP, SEXP posSXP, SEXP a1_SXP, SEXP a2_SXP, SEXP samples_SXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type get_info(get_infoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pos(posSXP);
    Rcpp::traits::input_parameter< CharacterVector >::type a1(a1_SXP);
    Rcpp::traits::input_parameter< CharacterVector >::type a2(a2_SXP);
    Rcpp::traits::input_parameter< CharacterVector >::type samples(samples_SXP);
    snp_filter F(chr, pos, a1, a2);
    rcpp_result_gen = Rcpp::wrap(read_vcf_filtered(filename, get_info, F, samples));
    return rcpp_result_gen;
END_RCPP
}
