#include <Rcpp.h>
#include <string>
// #include "gaston/flip_strand.h"
#include "flip_strand.h"

using namespace Rcpp;


//[[Rcpp::export]]
List check_alleles(CharacterVector A1, CharacterVector A2, CharacterVector B1, CharacterVector B2) {
  int n = A1.length();
  if(n != A2.length() || n != B1.length() || n != B2.length())
    stop("Dimensions mismatch");

  std::string a1, a2;
  LogicalVector ok(n), swap(n), flip(n);
  for(int i = 0; i < n; i++) {
     if( !strcmp( CHAR(A1[i]), CHAR(B1[i]) ) && !strcmp( CHAR(A2[i]), CHAR(B2[i]) ) ) {
       ok[i] = true;
       swap[i] = false;
       flip[i] = false;
       continue;
     }
     if( !strcmp( CHAR(A1[i]), CHAR(B2[i]) ) && !strcmp( CHAR(A2[i]), CHAR(B1[i]) ) ) {
       ok[i] = true;
       swap[i] = true;
       flip[i] = false;
       continue;
     }
     a1 = flip_strand( CHAR(A1[i]) );
     a2 = flip_strand( CHAR(A2[i]) );
     if(a1 == CHAR(B1[i]) && a2 ==  CHAR(B2[i]) ) {
       ok[i] = true;
       swap[i] = false;
       flip[i] = true;
       continue;
     }
     if(a1 == CHAR(B2[i]) && a2 ==  CHAR(B1[i]) ) {
       ok[i] = true;
       swap[i] = true;
       flip[i] = true;
       continue;
     }
     ok[i] = false;
  }
  List L;
  L["OK"] = ok;
  L["swap"] = swap;
  L["flip"] = flip;
  return L;
}

/*CharacterVector flip_strand(CharacterVector Allele) {
  std::vector<std::string> R;
  for(auto x : Allele) {
    R.push_back( flip_strand( CHAR(x) ) );
  }
  return wrap(R);
}*/

RcppExport SEXP gg_check_alleles(SEXP A1SEXP, SEXP A2SEXP, SEXP B1SEXP, SEXP B2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type A2(A2SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type B2(B2SEXP);
    rcpp_result_gen = Rcpp::wrap(check_alleles(A1, A2, B1, B2));
    return rcpp_result_gen;
END_RCPP
}

