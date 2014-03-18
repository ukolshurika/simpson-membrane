#include "bound.h"

#include <cmath>
#include <iostream>

#include "simpson.h"
#include "utils.h"

const double k2Sqrt3 = 2/sqrt(3);

struct SFunctor{
  SFunctor(){};
  double operator()(double x) const{
    return sqrt(1+25*pow(x, 10));
  }
};

struct HFunctor{
  HFunctor(const Bound& bound): b(bound){};

  double operator()(int i) const{
    return b.B1(i)/b.B2(i);
  }

  Bound b;
};

DBound::Bound(const Membrane& m):m_(m){};

double F(double x){
  return 1-x*x*x*x*x*x;
}

double DBound::operator()(int i ) const{
  // std::cerr << "CONSTRAINED" << std::endl;
  // std::cerr << B1(x)/B2(x) << std::endl;
  return B1(i)/B2(i)*pow(k2Sqrt3*H(x)/(m_.q_*Rho(i)), m_.n_);
}

double DBound::Rho(int i) const{
  return sqrt((xc[i]-x0[i])*(xc[i]-x0[i]) + (yc[i]-F(x0[i]))*(yc[i]-F(x0[i])));
}

double DBound::RhodRho(int i) const{
  return sqrt((xc[i+1]-x0[i+1])*(xc[i+1]-x0[i]) + (yc[i +1]-F(x0[i+1]))*(yc[i+1]-F(x0[i+1])));
}

double DBound::S(int i) const{
  return Simpson::Integrate(x0[i], Matrix::RZero(), kSimpsonStep, SFunctor())+
   Simpson::Integrate(0, x1[i], kSimpsonStep, SFunctor());
}

double DBound::SdS(int i) const{
  return S[i+1];
}

double DBound::Alpha(int i) const{
  return 2*asin( 0.5*sqrt((x1[i]-x0[i])*(x1[i]-x0[i]) + (F(x1[i])-F(x0[i]))*(F(x1[i])-F(x0[i])))/Rho[i]);
}

double DBound::AlphadAlpha(int i) const{
  return Alpha[i+1];
}

double DBound::B1(int i) const{
  // std::cerr << x << ' '<< Rho(x)*dAlpha(x) << " alpha:" << Alpha(x)<< ' ' <<dRho(x) << ' ' << dS(x) << std::endl;
  return Rho(i)*Alpha(i)+S(i)-B2(i);//Rho(x)*dAlpha(x)+Alpha(x)*dRho(x)+dS(x);
}

double DBound::B2(int i) const{
  return Rho(i)*Alpha(i)+S(i);
}

double DBound::SigmaE(int i) const{
  // double h = H(x);
  // std::cerr << m_.q_ <<  ' ' << Rho(x) << ' ' << h << std::endl;
  return sqrt(3)/2*m_.q_*Rho(i)/H(i)/0.02;
}

double DBound::H(int i) const{
  // std::cerr << Simpson::Integrate(x, Matrix::RZero(), kSimpsonStep, HFunctor(*this)) << std::endl;
  return m_.h1_*exp(-Simpson::Integrate2(0, i, kSimpsonStep, HFunctor(*this)));
}


