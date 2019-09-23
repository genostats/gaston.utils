#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "ms_files.h"

using namespace Rcpp;

msFile::msFile(std::string file, int rep) : filename(file), in( (char *) &filename[0u] ), rep(rep) {
  start();
}

msFile::msFile(const char * file, int rep) : filename(file), in(file), rep(rep) {
  start();
}

msFile::msFile(const CharacterVector Filename, int rep) : filename(Filename[0]), in( (const char *) &filename[0u] ), rep(rep) {
  start();
}

msFile::~msFile() {
  in.close();
}


void msFile::start() {
  if(!in.good())
    stop("Can't open file");

  std::string line;

  std::getline(in, line); // ligne de commande
  std::istringstream li(line);
  std::string s;
  if(!(li >> s >> nsamples))
    stop("ms file format error");


  for(int r = 0; r < rep; r++) {
    while(std::getline(in, line)) 
      if(line.substr(0,2) == "//")  
        break;
  }

  if(line.substr(0,2) != "//")
    stop("ms file format error");

  std::getline(in, line); // ligne segsites;
  li = std::istringstream(line);
  if(!(li >> s >> segsites))
    stop("ms file format error");

  std::getline(in, line); // ligne positions
  li = std::istringstream(line);
  if(!(li >> s))
    stop("ms file format error");
  double p;
  while(li >> p) 
    pos.push_back(p);
  
}

bool msFile::read_haplotype(std::vector<char> & haplo) {
  std::string line;
  if(!std::getline(in, line))
    return false;

  if(line.length() == 0)
    return false;

  if(line.length() != segsites)
    stop("ms file format error");

  haplo.resize(segsites);
  for(int i = 0; i < segsites; i++) 
    haplo[i] = line[i] - '0';

  return true;
};
  
bool msFile::read_genotype(std::vector<char> & geno) {
  std::string line;
  if(!std::getline(in, line))
    return false;

  if(line.length() == 0)
    return false;

  if(line.length() != segsites)
    stop("ms file format error");

  geno.resize(segsites);
  for(int i = 0; i < segsites; i++) 
    geno[i] = line[i] - '0';

  // deuxième haplotype
  if(!std::getline(in, line))
    return false;

  if(line.length() == 0)
    return false;

  if(line.length() != segsites)
    stop("ms file format error");

  for(int i = 0; i < segsites; i++) 
    geno[i] += line[i] - '0';

  return true;
}
