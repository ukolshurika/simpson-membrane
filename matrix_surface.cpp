#include "matrix_surface.h"

#include <cmath>

double MatrixSurface::operator () (double x) const {
  return -1*x*x+1;
}

double MatrixSurface::derivative (double x) const {
  return (MatrixSurface(x)+MatrixSurface(x+DELTA))/DELTA;
  // return -2*x;
}

// correct(???) for all matrixes
double MatrixSurface::alphaConstrained () const {
  return M_PI_2 - asinh( derivative(right_zero()));
}


//TODO  method to calculate zero automaticly ?!
double MatrixSurface::right_zero () const {
  return 1;
}

double MatrixSurface::normal(double x) const{
	return MatrixSurface(x)-1/(derivative(x))*(x-x);
}