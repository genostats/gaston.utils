// double emploi avec gaston m4_ld.cpp
// à rationnaliser un jour
#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "ld.h"

using namespace Rcpp;

double LD(matrix4 & A, double mu1, double mu2, double v, size_t x1, size_t x2) {
  double LD = 0;
  double gg[16];
  gg[3] = gg[7] = gg[11] = gg[12] = gg[13] = gg[14] = gg[15] = 0;
  gg[0] = (-mu1)*(-mu2);
  gg[1] = (-mu1)*(1.-mu2);
  gg[2] = (-mu1)*(2.-mu2);

  gg[4] = (1.-mu1)*(-mu2);
  gg[5] = (1.-mu1)*(1.-mu2);
  gg[6] = (1.-mu1)*(2.-mu2);

  gg[8] = (2.-mu1)*(-mu2);
  gg[9] = (2.-mu1)*(1.-mu2);
  gg[10]= (2.-mu1)*(2.-mu2);

  for(size_t i = 0; i < A.true_ncol; i++) {
    uint8_t g1 = A.data[x1][i];
    uint8_t g2 = A.data[x2][i];
    for(int ss = 0; ss < 4; ss++) {
      LD += gg[ ((g1&3)*4) + (g2&3) ];
      g1 >>= 2;
      g2 >>= 2;
    }
  }
  return LD/(v*(A.ncol-1));
}

