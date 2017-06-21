#include <Rcpp.h>
#include <iostream>
#include <string>
#ifndef TOKEN
#define TOKEN
using namespace Rcpp;

int token_position(std::string s, std::string token);

template<typename scalar_t>
scalar_t sto(const std::string & x);

template<typename scalar_t>
scalar_t token_at_position(std::string s, int pos) {
  std::istringstream ss(s);
  std::string token;
  for(int i = 0; i < pos && std::getline(ss, token, ':'); i++) {}
  std::getline(ss, token, ':');
  scalar_t r = sto<scalar_t>(token);
  return r;
}
#endif
