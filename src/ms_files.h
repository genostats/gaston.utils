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
  int nreplicates; // le nombre de replicats (tous de nsamples)

  std::vector<int> SEGSITES; // le nombre de Segsites, réplicat par réplicat
  int total_segsites; // nombre total de segsites (sommé sur tous les réplicats)

  int rep; // le replicat en cours de lecture
  int segsites; // nombre de SNPs du replicat en cours de lecture
  int snp; // numero du SNP du réplicat en cours de lecture... (de 0 à segsites-1)

  // will contain all data of the current replicate
  std::vector< std::vector<char> > currentRep;


  msFile(std::string file);
  msFile(const char * file);
  msFile(const CharacterVector Filename);
  ~msFile();


private:  
  void start(); // appelé par les constructeurs pour tout mettre en route
  bool next_replicate(); // lecture du réplicat suivant dans currentRep. False en fin de fichier.
public:
  // read one SNP for all individuals. All replicates are considered one after another.
  // sends true while there are some still SNPs to read.
  bool read_SNP(std::vector<char> & haplo);
};

#endif


