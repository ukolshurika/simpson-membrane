#include "ideal_sliding.h"

#include <cmath>
#include <iostream>

#include "simpson.h"

int kSimpsonStep = 999;

double Circle::center(double a, double alpha){
  return -1*a/tan(alpha);
}

IdealSliding::IdealSliding(const MatrixSurface& ms, const Membrane& m, double h1): ms_(ms), m_(m), h1_(h1){
}

double IdealSliding::operator () (double x) const {
  return B1(x)/B2(x)*pow((2*h(x)*m_.sigma_b_)/(sqrt(3)*m_.q_*Rho(x) - 2), m_.n_);
}

double IdealSliding::Alpha(double x) const{
  return M_PI_2 - atan(ms_.dNormal(x));
}

double IdealSliding::dAlpha(double x) const{
  return (Alpha(x+DELTA) - Alpha(x))/DELTA;
}

double IdealSliding::S(double x) const{
  //WARN: integration order is changed!
  return Simpson::Integrate(x, ms_.RightZero(), kSimpsonStep, SFunctor((*this)));
}

double IdealSliding::dS(double x) const{
  return (S(x+DELTA) - S(x))/DELTA;
}

double IdealSliding::Rho(double x) const{
  return (ms_(x) - Circle::center(m_.a_, Alpha(x)))*(ms_(x) - Circle::center(m_.a_, Alpha(x)))+x*x;
  // return sqrt((ms_(x) - Circle::center(m_.a_, Alpha(x)))*(ms_(x) - Circle::center(m_.a_, Alpha(x)))+x*x);
  // return 0.1;
}

double IdealSliding::dRho(double x) const{
  return (Rho(x+DELTA) - Rho(x))/DELTA;
}

double IdealSliding::B1(double x) const {
  return Rho(x) * dAlpha(x) + Alpha(x)*dRho(x) + dS(x);
}

double IdealSliding::B2(double x) const{
  return Rho(x)*Alpha(x) + S(x);
}

double IdealSliding::h(double x) const{
  // return 0.1;
  return -1*h1_*exp(Simpson::Integrate(x, ms_.RightZero(), kSimpsonStep, HFunctor((*this))));
}

HFunctor::HFunctor (const IdealSliding& is): is_(is){};
SFunctor::SFunctor (const IdealSliding& is): is_(is){};

double HFunctor::operator () (double x) const{
  return is_.B1(x)/is_.B2(x);
}

double SFunctor::operator () (double x) const{
  // return 0.1;
  return sqrt(1+is_.ms_.SecondDerivative(x)*is_.ms_.SecondDerivative(x));
}