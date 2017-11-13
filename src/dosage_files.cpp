#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "dosage_files.h"
#include "gzstream.h"
#include "token.h"
#include "read_gen_line.h"
#include "read_vcf_line.h"

using namespace Rcpp;

dosages::dosages(std::string file) : filename(file), in( (char *) &filename[0u] ) {
  start();
}

dosages::dosages(const char * file) : filename(file), in( file ) {
  start();
}

dosages::dosages(const CharacterVector Filename) : filename(Filename[0]), in( (const char *) &filename[0u] ) {
  start();
}

dosages::~dosages() {
  in.close();
}

// se met en début de fichier 
// initialise le vectuer samples quand c'est un VCF
void dosages::start() {
  if(!in.good()) 
    stop("Can't open file");
  if(!std::getline(in, line)) 
    stop("File is empty");

  if(line.substr(0,1) == "#") { // VCF
    // Rcout << "VCF\n";
    type = 10;
    // skip description informations
    while(std::getline(in, line)) {
    if(line.substr(0,1) != "#") stop("Bad VCF format");
      if(line.substr(0,2) != "##") {
        read_vcf_samples(line, samples);
        break; // fin
      }
    }
    // read first data line
    if(std::getline(in, line)) 
      good = true;
    else
      good = false;

  } else if(line.substr(0,3) == "---") { // Impute2
    // Rcout << "Impute2\n";
    type = 1;
    good = true;
  } else {
    in.close();
    stop("Unknown file format");
  }
}

bool dosages::read_line(std::vector<double> & dosage, std::string & snp_id, 
                   int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  if(!good) return false;
  if(type == 1) {
    chr = "NA";
    parse_gen_line<double>(line, dosage, snp_id, snp_pos, A1, A2);
  }
  if(type == 10)
    parse_vcf_line_dosages<double>(line, dosage, snp_id, snp_pos, chr, A1, A2);
  // read next line now...
  if(std::getline(in, line))
    good = true;
  else
    good = false;

  return true;
}

bool dosages::read_line(std::vector<float> & dosage, std::string & snp_id, 
                   int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  if(!good) return false;
  if(type == 1) {
    chr = "NA";
    parse_gen_line<float>(line, dosage, snp_id, snp_pos, A1, A2);
  }
  if(type == 10)
    parse_vcf_line_dosages<float>(line, dosage, snp_id, snp_pos, chr, A1, A2);
  // read next line now...
  if(std::getline(in, line))
    good = true;
  else
    good = false;

  return true;
}


