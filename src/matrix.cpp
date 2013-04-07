#include "matrix.h"

#include <cassert>
#include <cmath>

const double Matrix::kK = 1.5;
const double Matrix::kB = 4.5;
const double Matrix::kKBSquare = kK*kK*kB*kB;

double Matrix::operator()(double x){
  return kB*(1-pow(x, kK));
  //return ;//-1.1918*x+1.1918;
}

// static
double Matrix::Derivative(double x){
  return -kK*pow(x, kK-1);//-1.1918;
}

// static
double Matrix::RZero(){
  return 1.0;
}