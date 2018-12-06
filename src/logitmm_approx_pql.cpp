#include <Rcpp.h>
#include "approx_pql.h"
using namespace Rcpp;


//[[Rcpp::export]]
List logitmm_approxpql(NumericVector Y, NumericMatrix X, NumericVector Sigma, NumericMatrix U, double tol, bool verbose) {

  approx_pql<double> A(as<Eigen::VectorXd>(Y), as<Eigen::MatrixXd>(X), as<Eigen::VectorXd>(Sigma), as<Eigen::MatrixXd>(U), tol);
  A.optimize(0, verbose); 

  List L;
  L["h2"] = A.h2;
  L["v"] = A.dlmm_object.v;
  L["beta"] = A.Beta;
  L["varbeta"] = A.VarBeta;
  return L;
}


RcppExport SEXP gg_logitmm_approxpql(SEXP YSEXP, SEXP XSEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP tolSEXP, SEXP verb) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verb);
    rcpp_result_gen = Rcpp::wrap(logitmm_approxpql(Y, X, Sigma, U, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}

