#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "token.h"
#include "read_gen_line.h"
#include "read_vcf_line.h"

using namespace Rcpp;

class dosages {
public:
  std::string filename;
  igzstream in;
  std::string line;
  int type;
  bool good;
  
  dosages(std::string file);
  dosages(const char * file);
  dosages(const CharacterVector Filename);
  ~dosages();


private:  
  void start();

public:
  bool read_line(std::vector<double> & dosage, std::string & snp_id, 
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

  bool read_line(std::vector<float> & dosage, std::string & snp_id, 
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

};

