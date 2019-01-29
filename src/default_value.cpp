#include "default_value.h"
#include <string>
#include <Rcpp.h>

template<>
double default_value<double>() {
  return NA_REAL;
}

template<>
float default_value<float>() {
  return NA_REAL;
}

template<>
int default_value<int>() {
  return NA_INTEGER;
}

template<>
std::string default_value<std::string>() {
  return "";
}

