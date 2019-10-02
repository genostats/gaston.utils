#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "ms_files.h"

using namespace Rcpp;

msFile::msFile(std::string file) : filename(file), in( (char *) &filename[0u] ) {
  start();
}

msFile::msFile(const char * file) : filename(file), in(file) {
  start();
}

msFile::msFile(const CharacterVector Filename) : filename(Filename[0]), in( (const char *) &filename[0u] ) {
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
  if(!(li >> s >> nsamples >> nreplicates))
    stop("ms file format error");

  SEGSITES.reserve(nreplicates);
  total_segsites = 0;
  // premiere passe sur le fichier.
  // on lit tous les réplicats...
  for(int r = 0; r < nreplicates; r++) {
    // aller au début d'un replicat
    while(std::getline(in, line)) 
      if(line.substr(0,2) == "//")  
        break;
    // vérifier que tout va bien
    if(line.substr(0,2) != "//")
      stop("ms file format error");

    // lire ligne segsites
    std::getline(in, line);
    li = std::istringstream(line);
    if(!(li >> s >> segsites))
      stop("ms file format error");

    SEGSITES.push_back(segsites);
    total_segsites += segsites;
  }

  // rewind file
  in.seekg(0);

  // Go to first replicate
  rep = 0;
  if(!next_replicate())
    stop("ms file format error"); 
}

bool msFile::next_replicate() {
  std::string line;
  std::istringstream li(line);
  std::string s;
  // aller au début d'un replicat
  while(std::getline(in, line)) 
    if(line.substr(0,2) == "//")
      break;

  // vérifier qu'on y est bien !
  if(line.substr(0,2) != "//")
    return false;

  // lire lignes segsites et positions 
  // on ne teste pas la présence possible d'un arbre (option -T) ou d'une proba (-t et -s à la fois)
  std::getline(in, line);
  li = std::istringstream(line);
  if(!(li >> s >> segsites))
    stop("ms file format error");
  // ligne position
  std::getline(in, line);
  // On lit tout 
  currentRep.clear();
  while( std::getline(in, line) ) {
    if(line.length() == 0) // fin du réplicat
      break;
    if(line.length() != segsites)
      stop("ms file format error");
    std::vector<char> haplo(segsites);
    for(int i = 0; i < segsites; i++) 
      haplo[i] = line[i] - '0';
    currentRep.push_back(haplo);
  }
  rep++;
  snp = 0;
  return true;
}


bool msFile::read_SNP(std::vector<char> & v) {
  v.clear();
  if(snp == segsites) // fin du replicat en cours
    if(!next_replicate())
      return false;
  
  for(int i = 0; i < nsamples; i++) {
    v.push_back( currentRep[i][snp] );
  }
  snp++;
  return true;
};
 
