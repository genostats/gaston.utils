#include <Rcpp.h>
#include <iostream>
#include <fstream>
#ifndef msFiles
#define msFiles

using namespace Rcpp;

class msFile {
public:
  std::string filename;
  std::ifstream in;
  int nsamples; // le nombre de 'samples' (chromosomes) générés par ms
  int rep; // le replicat à lire
  bool good;
  int segsites; // nombre de SNPs
  std::vector<double> pos; // la ligne 'positions'

  msFile(std::string file, int rep = 1);
  msFile(const char * file, int rep = 1);
  msFile(const CharacterVector Filename, int rep = 1);
  ~msFile();


private:  
  void start();

public:
  bool read_haplotype(std::vector<char> & haplo);
  bool read_genotype(std::vector<char> & geno);

};

#endif


