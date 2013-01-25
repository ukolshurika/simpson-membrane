#include "matrix_surface.h"

#include <cmath>

double MatrixSurface::operator () (double x) const {
  return -1*x*x+1;
}

double MatrixSurface::Derivative (double x) const {
  return ((*this)(x)+(*this)(x+DELTA))/DELTA;
  // return -2*x;
}

double MatrixSurface::SecondDerivative (double x) const {
  return ((*this)(x+DELTA) - 2*(*this)(x) + (*this)(x-DELTA))/DELTA/DELTA;
  // return -2;
}

// correct(???) for all matrixes
double MatrixSurface::AlphaConstrained () const {
  return M_PI_2 - asinh( Derivative(RightZero()));
}


//TODO  method to calculate zero automaticly ?!
double MatrixSurface::RightZero () const {
  return 1;
}

double MatrixSurface::Normal(double x) const{
	return (*this)(x)-1/(Derivative(x))*(x-x);
}

double MatrixSurface::dNormal(double x) const{
	return (*this)(x)-1/(Derivative(x))*(x-x);
}