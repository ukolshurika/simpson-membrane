#include "free_deformation.h"

#include <cmath>

FreeDeformation::FreeDeformation(double h0, double q, double n): h0_(h0), q_(q), n_(n) {}

double FreeDeformation::operator () (double alpha) const {
  return (1/alpha-1/tan(alpha))*pow((2*h0_*sin(alpha)*sin(alpha)/(sqrt(3)*q_*alpha) -1), n_);
}

double FreeDeformation::h(double alpha){
  return sin(alpha)/alpha;
};

