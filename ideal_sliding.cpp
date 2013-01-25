#include "ideal_sliding.h"

#include <cmath>

#include "simpson.h"

#define SIMSON_STEP 999 

double Cyrcle::center(double a, double alpha){
  return -1*a/tan(alpha);
}

IdealSliding::IdealSliding(const MatrixSurface& ms, const Membrane& m, double h1): ms_(ms), m_(m), h1_(h1){
}

double IdealSliding::operator () (double x) const {
  return B1(x)/B2(x)*pow((2*h(x)*m_.sigma_b_)/(sqrt(3)*m_.q_*Rho(x) - 2), m_.n_);
}

double IdealSliding::Alpha(double x) const{
  return M_PI_2 - atanh(ms_.dNormal(x));
}

double IdealSliding::dAlpha(double x) const{
  return (Alpha(x) + Alpha(x+DELTA))/DELTA;
}

double IdealSliding::S(double x) const{
  //WARN: integration order is changed!
  return -1*Simpson::Integrate(1, x, SIMSON_STEP, SFunctor((*this)));
}

double IdealSliding::dS(double x) const{
  return (S(x) + S(x+DELTA))/DELTA;
}

double IdealSliding::Rho(double x) const{
  return sqrt((m_(x) - Cyrcle::center(m_.a_, Alpha(x) )) * (m_(x) - Cyrcle::center(m_.a_, Alpha(x)))+x*x);
}

double IdealSliding::dRho(double x) const{
  return (Rho(x) + Rho(x+DELTA))/DELTA;
}

double IdealSliding::B1(double x) const {
  return Rho(x) * dAlpha(x) + Alpha(x)*dRho(x) + dS(x);
}

double IdealSliding::B2(double x) const{
  return Rho(x)*Alpha(x) + S(x);
}

double IdealSliding::h(double x) const{
  return h1_*exp(Simpson::Integrate(1, x, SIMSON_STEP, HFunctor((*this))));
}

HFunctor::HFunctor (const IdealSliding& is): is_(is){};
SFunctor::SFunctor (const IdealSliding& is): is_(is){};

double HFunctor::operator () (double x) const{
  return is_.B1(x)/is_.B2(x);
}

double SFunctor::operator () (double x) const{
  return sqrt(1+is_.ms_.SecondDerivative(x)*is_.ms_.SecondDerivative(x));
}