#include "ideal_sliding.h"

#include <cmath>
#include <iostream>
#include <cassert>

#include "simpson.h"
#include "utils.h"


namespace {
const double kEpsilon = 1e-9;
const int kSimpsonStep = 999;
bool eql(double u, double v){
  return abs(u - v) < kEpsilon;
}

}

double Circle::center(double a, double alpha){
  return -1*a/tan(alpha);
}

IdealSliding::IdealSliding(const MatrixSurface& ms, const Membrane& m): ms_(ms), m_(m){
  h1_=Alpha(ms_.RightZero());
}

double IdealSliding::operator () (double x) const {
  double y = B1(x)/B2(x)*pow((2*h(x)*m_.sigma_b_)/(sqrt(3)*m_.q_*Rho(x)) - 2, m_.n_);
  // DCHECK(B1(x)/B2(x)>=0);
  return y;
}

double IdealSliding::Alpha(double x) const{
  double y = M_PI_2 - atan(ms_.dNormal(x));
  DCHECK(y>=0);
  return y;
}

double IdealSliding::dAlpha(double x) const{
  double y = (Alpha(x+DELTA) - Alpha(x))/DELTA;
  // assert(y>=0);
  return y;
}

double IdealSliding::S(double x) const{
  //WARN: integration order is changed!
  double y = Simpson::Integrate(x, ms_.RightZero(), kSimpsonStep, SFunctor((*this)));
  if(eql(x, ms_.RightZero()))
    y=0;
  DCHECK(y>=0);
  assert(y>=0);
  return y;
}

double IdealSliding::dS(double x) const{
  double y = SFunctor((*this))(x);
  // double y = (S(x+DELTA) - S(x))/DELTA;
  DCHECK(y>=0);
  return y;
}

double IdealSliding::Rho(double x) const{
  double y = sqrt((ms_(x) - Circle::center(x, Alpha(x)))*(ms_(x) - Circle::center(x, Alpha(x)))+x*x);
  DCHECK(y>=0);
  return y;
  // return 0.1;
}

double IdealSliding::dRho(double x) const{
  double y = (Rho(x+DELTA) - Rho(x))/DELTA;
  // DCHECK(y>=0);
  return y; 
}

double IdealSliding::B1(double x) const {
  double y = Rho(x) * dAlpha(x) + Alpha(x)*dRho(x) + dS(x);
  // std::cout << x << ": "<< Rho(x) * dAlpha(x) << ' ' << Alpha(x)*dRho(x) << ' ' << dS(x) << std::endl;
  // DCHECK(y>=0);
  return y;
}

double IdealSliding::B2(double x) const{
  double y = Rho(x)*Alpha(x) + S(x);
  DCHECK(y>=0);
  return y;
}

double IdealSliding::h(double x) const{
  // return 0.1;
  // double y = h1_*exp(-Simpson::Integrate(x, ms_.RightZero(), kSimpsonStep, HFunctor((*this))));
  double y = HFunctor(*this)(x);
  // std::cout << x <<' ' <<y<< std::endl;
  // DCHECK(y>=0);
  return y;
}

double IdealSliding::SigmaE(double x) const{
  return sqrt(3)/2*m_.q_*Rho(x)/h(x);
}

HFunctor::HFunctor (const IdealSliding& is): is_(is){};
SFunctor::SFunctor (const IdealSliding& is): is_(is){};

double HFunctor::operator () (double x) const{
  // std::cout << x << ": "<< is_.B1(x) << " " << is_.B2(x) << std::endl;
  return is_.B1(x)/is_.B2(x);

}

double SFunctor::operator () (double x) const{
  // return 1;
  return sqrt(1+is_.ms_.Derivative(x)*is_.ms_.Derivative(x));
}