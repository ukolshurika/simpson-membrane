#include "matrix_surface.h"

#include <cmath>

namespace{
  const double kDelta = 1e-3;
  const double kInvDelta = 1/kDelta;
  const double kInvSquareDelta = kInvDelta*kInvDelta;
}

double MatrixSurface::operator () (double x) const {
  return -3*x+3/*-pow(x, 1.5)+1*/;
}

double MatrixSurface::Derivative (double x) const {
  // return ((*this)(x+kDelta) - (*this)(x))*kInvDelta;
  return -3/*-1.5*sqrt(x)*/;
}

double MatrixSurface::SecondDerivative (double x) const {
  // return ((*this)(x+kDelta) - 2*(*this)(x) + (*this)(x-kDelta))*kInvSquareDelta;
  return 0/*-1.5/sqrt(x)*/;
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