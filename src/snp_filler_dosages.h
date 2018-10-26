#include <Rcpp.h>
#include "gzstream.h"
#include "dosage_files.h"
#include "snp_filler.h"

#ifndef GASTON_SNP_FILLER_DOSAGES
#define GASTON_SNP_FILLER_DOSAGES

using namespace Rcpp;

template<typename scalar_t>
class snp_filler_dosages : public snp_filler<scalar_t> {
  public:

  dosages in;

  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> SNP_ID, CHR, AL1, AL2;
  std::vector<int> POS;
  std::vector<scalar_t> dosage;
  std::vector<double> F1;
  std::vector<double> F2;

  int beg, end;
  int nb_snps;
  int i;

  snp_filler_dosages(CharacterVector filename, int beg_, int end_, int n)
    : snp_filler<scalar_t>(), in(filename), beg(beg_), end(end_), nb_snps(n), i(0) {
  };

  bool snp_fill(scalar_t * SNP) {
    this->monomorphic = true;
    if(!in.read_line(dosage, snp_id, snp_pos, chr, A1, A2)) {
      dosage.clear();
      return false; // EOF
    }

    i++;
    if(i < beg) {
      dosage.clear();
      return true;
    }
    if(i > end) {
      dosage.clear();
      return false; 
    }
    if(dosage.size() != nb_snps) {
      dosage.clear()
      return false;
    }

    SNP_ID.push_back(snp_id);
    POS.push_back(snp_pos);
    CHR.push_back(chr);
    AL1.push_back(A1);
    AL2.push_back(A2);

    // ajout deux cols fréquences et remplissage SNP
    // double ss = std::accumulate(dosage.begin(), dosage.end(), 0.0)/dosage.size()/2.0;
    scalar_t ss = 0;
    for(int i = 0; i < nb_snps; i++) {
       ss += SNP[i] = dosage[i];
    }
    ss /= 2.0;

    F1.push_back(1.0 - ss);
    F2.push_back(ss);

    if(ss < 1.0 || ss > 0.0) {
       this->monorphic = false;
    }
    return true;
  }
};

