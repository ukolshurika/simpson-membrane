#include "bound.h"

#include <cmath>
#include <iostream>

#include "simpson.h"
#include "utils.h"

const double k2Sqrt3 = 2/sqrt(3);

struct SFunctor{
  SFunctor(){};
  double operator()(double x) const{
    return sqrt(1+Matrix::kKBSquare*x*x);
  }
};

struct HFunctor{
  HFunctor(const Bound& bound): b(bound){};

  double operator()(double x) const{
    return b.B1(x)/b.B2(x);
  }

  Bound b;
};

Bound::Bound(const Membrane& m):m_(m){};

double Bound::operator()(double x) const{
  // std::cerr << "CONSTRAINED" << std::endl;
  // std::cerr << B1(x)/B2(x) << std::endl;
  return B1(x)/B2(x)*pow(k2Sqrt3*H(x)/(m_.q_*Rho(x)), m_.n_);
}

double Bound::Rho(double x) const{
  return sqrt(x*x +x*x/(Matrix::kKBSquare*pow(x, 2*Matrix::kK-2)));
}

double Bound::dRho(double x) const{
  return (x*x+(2-Matrix::kK)*pow(x, 3-2*Matrix::kK)/Matrix::kKBSquare)/Rho(x);
}

double Bound::S(double x) const{
  return Simpson::Integrate(x, Matrix::RZero(), kSimpsonStep, SFunctor());
}

double Bound::dS(double x) const{
  return -sqrt(1+Matrix::kKBSquare*pow(x, 2*Matrix::kK-2));
}

double Bound::Alpha(double x) const{
  return M_PI_2 - atan(1/(Matrix::kK*Matrix::kB*pow(x, Matrix::kK-1)));
}
double Bound::dAlpha(double x) const{
  double k = Matrix::kK;
  double b = Matrix::kB;
  return (k*(k-1)*b*pow(x, k-2)/(1+Matrix::kKBSquare*pow(x, 2*k-2)));
}

double Bound::B1(double x) const{
  // std::cerr << x << ' '<< Rho(x)*dAlpha(x) << " alpha:" << Alpha(x)<< ' ' <<dRho(x) << ' ' << dS(x) << std::endl;
  return Rho(x)*dAlpha(x)+Alpha(x)*dRho(x)+dS(x);
}

double Bound::B2(double x) const{
  return Rho(x)*Alpha(x);
}

double Bound::SigmaE(double x) const{
  // double h = H(x);
  // std::cerr << m_.q_ <<  ' ' << Rho(x) << ' ' << h << std::endl;
  return sqrt(3)/2*m_.q_*Rho(x)/H(x)/0.02;
}

double Bound::H(double x) const{
  // std::cerr << Simpson::Integrate(x, Matrix::RZero(), kSimpsonStep, HFunctor(*this)) << std::endl;
  return m_.h1_*exp(Simpson::Integrate(x, Matrix::RZero(), kSimpsonStep, HFunctor(*this)));
}


