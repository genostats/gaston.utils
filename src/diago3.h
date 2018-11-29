#include <RcppEigen.h>
#include <iostream>
#include "gaston/matrix-varia.h"
#include "gaston/optimize_2.h"
#ifndef GASTONdiag_lmm_likelihood
#define GASTONdiag_lmm_likelihood

using namespace Rcpp;
using namespace Eigen;

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;


// structure qui hérite fun pour l'optimisation par Brent
// Normalement l'utilisation doit être
// MATRIX = MatrixXd, VECTOR = VectorXd, scalar_t = double
// MATRIX = MatrixXf, VECTOR = VectorXf, scalar_t = float
template<typename scalar_t>
class diag_lmm_likelihood : public fun<scalar_t> {
  public:
    int p, n, r;
    const VECTOR<scalar_t> Y;
    MATRIX<scalar_t> X;
    const MATRIX<scalar_t> Sigma;
    VECTOR<scalar_t> P0y;
    scalar_t v;
    MATRIX<scalar_t> XViX_i;
    VECTOR<scalar_t> Deltab;
    scalar_t d, log_d;
    VECTOR<scalar_t> V0b, V0bi;
    MATRIX<scalar_t> ViX, XViX, xtx;
    scalar_t yP0y;
    diag_lmm_likelihood(int p_, const VECTOR<scalar_t> & Y_, const MATRIX<scalar_t> & X_, const VECTOR<scalar_t> & Sigma_) : 
       p(p_), n(Sigma_.rows()), r(X_.cols()), Y(Y_), X(X_), Sigma(Sigma_) {
         Deltab = Sigma.bottomRows(n-p) - VECTOR<scalar_t>::Ones(n-p);
         XViX_i = MATRIX<scalar_t>(r,r);
    } ;

    void update(scalar_t h2) {
      V0b = h2*Sigma.bottomRows(n-p) + (1-h2)*VECTOR<scalar_t>::Ones(n-p); 
      V0bi = V0b.cwiseInverse();
      ViX = V0bi.asDiagonal() * X.bottomRows(n-p);  //      V_{0b}^{-1} X_b
      XViX = X.bottomRows(n-p).transpose() * ViX;   // X_b' V_{0b}^{-1} X_b
      sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

      P0y = V0bi.asDiagonal() * Y.tail(n-p) - ViX * ( XViX_i * ( ViX.transpose() * Y.tail(n-p) ) ) ;
      yP0y = P0y.dot( Y.tail(n-p).col(0) );
      v = yP0y / (n-r-p);
    }

    scalar_t f(scalar_t h2) {
      update(h2);
      return -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y) + (n-r-p)*(1-log((scalar_t)(n-r-p))) ); // avec le terme constant
    }

    // vraisemblance au dernier point calculé
    scalar_t likelihood() {
      return -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y) + (n-r-p)*(1-log((scalar_t)(n-r-p))) ); // avec le terme constant
    }

    scalar_t likelihood(scalar_t tau, scalar_t s2) {
      scalar_t vv = s2 + tau;
      scalar_t h2 = tau/vv;
      update(h2);
      return -0.5*(V0b.array().log().sum() + log_d + yP0y/vv + (n-r-p)*log(vv));
    }
 
    void df_ddf(scalar_t h2, scalar_t & df, scalar_t & ddf) {
      update(h2);
      VECTOR<scalar_t> DeltaP0y = Deltab.cwiseProduct( P0y );
      scalar_t y_PDP_y = P0y.dot( DeltaP0y );
  
      // VECTOR P0DeltaP0y = P0b * DeltaP0y;
      VECTOR<scalar_t> P0DeltaP0y  = V0bi.asDiagonal() * DeltaP0y - ViX * ( XViX_i * ( ViX.transpose() * DeltaP0y ) ) ;
      scalar_t y_PDPDP_y = DeltaP0y.dot( P0DeltaP0y );

      // scalar_t f = -0.5*(V0b.array().log().sum() + log_d + (n-r-p)*log(yP0y)+ (n-r-p)*(1-log((scalar_t)(n-r-p))));
      // **** calcul dérivées
      // MATRIX P0b = V0bi.asDiagonal();
      // P0b.noalias() -= ViX * XViX_i * ViX.transpose();
      // MATRIX DeltaP0 = Deltab.asDiagonal() * P0b;
      // scalar_t df0 = DeltaP0.trace() - (n-r-p)*y_PDP_y/yP0y;
      // scalar_t ddf0 = -trace_of_product(DeltaP0, DeltaP0) + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y));

      // alternative qui évite le calcul de P0b et Delta P0b (qui sont rapides à calculer
      // mais de grande taille ->O(n²r) remplacé par O(r²n) si je ne m'abuse
      MATRIX<scalar_t> B = ViX.transpose() * Deltab.asDiagonal() * ViX;
      MATRIX<scalar_t> C = ViX.transpose() * (Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab))).asDiagonal() * ViX;
      MATRIX<scalar_t> AB = XViX_i * B;
      scalar_t trace_PD = Deltab.cwiseProduct(V0bi).sum() - AB.trace();
      scalar_t trace_PDPD = Deltab.cwiseProduct(V0bi.cwiseProduct(Deltab.cwiseProduct(V0bi))).sum() - 2*trace_of_product(XViX_i, C) + trace_of_product(AB,AB);
  
      df = -0.5*(trace_PD - (n-r-p)*y_PDP_y/yP0y);
      ddf = -0.5*(-trace_PDPD + (n-r-p)*( 2*y_PDPDP_y/yP0y - (y_PDP_y*y_PDP_y)/(yP0y*yP0y)));
    }



    // *********** CALCUL DES BLUPS ************************
    // Attention P0y n'est que (P0y)b, les n-p dernières composantes ! (les p premières sont nulles)
    void blup(scalar_t h2, VECTOR<scalar_t> & beta, VECTOR<scalar_t> & omega, bool updateh2, bool only_first_r) {
      if(updateh2) update(h2);

      VECTOR<scalar_t> Sigmab = Sigma.bottomRows(n-p);
      omega = h2 * Sigmab.asDiagonal() * P0y;

      VECTOR<scalar_t> z = Y;
      z.tail(n-p) -= omega + (1-h2)*P0y;
      // Xb' Xb
      // MATRIX xtx( MATRIX(r,r).setZero().selfadjointView<Lower>().rankUpdate( X.bottomRows(n-p).transpose() ));
      xtx = MATRIX<scalar_t>(r,r).setZero();
      SelfAdjointView<MATRIX<scalar_t>, Lower> xtx_sa(xtx);
      xtx_sa.rankUpdate( X.bottomRows(n-p).transpose() );  
      MATRIX<scalar_t> xtx0( xtx_sa );
      MATRIX<scalar_t> xtxi(r,r); // et son inverse
      scalar_t d1, ld1;
      sym_inverse(xtx0, xtxi, d1, ld1, 1e-5); // détruit xtx0

      if(only_first_r) {
        beta = VECTOR<scalar_t>(r);
        beta = xtxi * X.bottomRows(n-p).transpose() * z.bottomRows(n-p);
      } else {
        beta = VECTOR<scalar_t>(r+p);
        beta.topRows(r) = xtxi * X.bottomRows(n-p).transpose() * z.bottomRows(n-p);
        beta.bottomRows(p) = z.topRows(p) - X.topRows(p) * beta.topRows(r);
     }
   }

};
#endif
