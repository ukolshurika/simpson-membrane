#include "matrix_surface.h"

#include <cmath>

namespace{
  const double kDelta = 1e-3;
  const double kInvDelta = 1/kDelta;
  const double kInvSquareDelta = kInvDelta*kInvDelta;
}

double MatrixSurface::operator () (double x) const {
  return -x*x+1;
}

double MatrixSurface::Derivative (double x) const {
  // return ((*this)(x+kDelta) - (*this)(x))*kInvDelta;
  return -2*x;
}

double MatrixSurface::SecondDerivative (double x) const {
  // return ((*this)(x+kDelta) - 2*(*this)(x) + (*this)(x-kDelta))*kInvSquareDelta;
  return -2;
}

// correct(???) for all matrixes
double MatrixSurface::AlphaConstrained () const {
  return M_PI_2 - asinh( Derivative(RightZero()));
}


//TODO  method to calculate zero automaticly ?!
double MatrixSurface::RightZero () const {
  return 1;
}


double MatrixSurface::dNormal(double x) const{
  return 1/(Derivative(x));
}