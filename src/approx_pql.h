#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston/logit_model.h"
#include "diago3.h"
#include <cmath>

using namespace Rcpp;
using namespace Eigen;

template<typename scalar_t>
scalar_t expit(scalar_t x) {
  return 1/(1+exp(-x));
}

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;


template<typename scalar_t>
class approx_pql {
  public:
  int n, r;  // nb individus, nb cols X

  VECTOR<scalar_t> Y;
  MATRIX<scalar_t> X; 
  VECTOR<scalar_t> Sigma;
  MATRIX<scalar_t> U; 
  double tol;

  VECTOR<scalar_t> Beta, Beta0;
  MATRIX<scalar_t> VarBeta; 
  diag_lmm_likelihood<scalar_t> dlmm_object;

  VECTOR<scalar_t> nu, pi_hat, Z, Omega;
  scalar_t h2;

  approx_pql(const VECTOR<scalar_t> Y_, const MATRIX<scalar_t> X_, const VECTOR<scalar_t> Sigma_, const MATRIX<scalar_t> U_, double tol) :
    n(Y_.rows()), r(X_.cols()), Y(Y_), X(X_), Sigma(Sigma_), U(U_), tol(tol),  
    Beta(r), Beta0(r), VarBeta(r,r), dlmm_object(0, Y, X, Sigma) {
    // Note : l'objet lmm a X et Y mal initialisés [sera updaté par U'X et U'Z dans optimize]
/*
    Rcout << "Init approx pql\n";
    Rcout << "Y = " ;
    for(int i = 0; i < 10 ; i++) Rcout << Y[i] << " ";
    Rcout << "...\n";
    Rcout << "X = \n";
    for(int i = 0; i < 6 ; i++) {
      for(int j = 0; j < X.cols(); j++) Rcout << X(i,j) << " ";
      Rcout << "\n";
    }
    Rcout << "...\n";
*/
  }

  void optimize(double h2_, bool verbose) {

/*
    Rcout << "Optimize approx pql\n";
    Rcout << "Y = " ;
    for(int i = 0; i < 10 ; i++) Rcout << Y[i] << " ";
    Rcout << "...\n";
    Rcout << "X = \n";
    for(int i = 0; i < 6 ; i++) {
      for(int j = 0; j < X.cols(); j++) Rcout << X(i,j) << " ";
      Rcout << "\n";
    }
    Rcout << "...\n";
*/

    h2 = h2_;
    // initialisation Beta par régression logistique classique
    logistic_model2<scalar_t>(Y, X, Beta0, VarBeta);
    Beta = Beta0;

    // Calcul de U.transpose * X  -> la matrice de covariable du LMM
    // le faire ici permet d'être OK avec le fait que X a pu changer 
    // entre la création de l'objet et l'appel d'optimize
    dlmm_object.X = U.transpose() * X;
      
    while(true) { 
      if(verbose) Rcout << "Beta = " << Beta.transpose() << " : ";
      // Calcul de Z = working vector  
      nu = X * Beta;
      pi_hat = nu.unaryExpr(std::ptr_fun(expit<scalar_t>));  // expit(nu), terme à terme
      Z = nu + (Y - pi_hat).cwiseQuotient( pi_hat.cwiseProduct( VECTOR<scalar_t>::Ones(n) - pi_hat ) ); // nu + (y - pi^) / ( pi^ (1-pi^) )
      // Calcul de Z1 = U.transpose * Z  -> le nouveau vecteur de réponse du LMM
      dlmm_object.Y = U.transpose() * Z;
      // A single step of Newton optimisation ...
      bool converged = dlmm_object.newton_max(h2, 0, 0.99, tol, 1, verbose, false);
      // compute blups for current value of h2
      dlmm_object.blup(h2, Beta, Omega, false, true);
      if(converged)
        break;
    }
    if(h2 == 0)  // fallback to logistic regression
      Beta = Beta0;  // Pas touche à VarBeta calculé par logistic_model
    else
      VarBeta = dlmm_object.v * dlmm_object.XViX_i;

  }

};
