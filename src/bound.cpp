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
    if(b.ordinate_=='y')
      return b.B1(x)/b.B2(x);
    if(b.ordinate_=='x')
      return (2-M_PI_2)/(Bound::kB+(M_PI_2-1))+(2-M_PI_2)*x;
  }

  Bound b;
};

const double Bound::kB = 4.5;

Bound::Bound(const Membrane& m, char ordinate):m_(m), ordinate_(ordinate){};

double Bound::operator()(double x) const{
  // std::cerr << "CONSTRAINED" << std::endl;
  // std::cerr << "!!!!"<< B1(x) << ' ' << B2(x) << std::endl;
  double multiplier;
  multiplier = (ordinate_ == 'y') ? 1/(x+M_PI_2/2) : ((2-M_PI_2)/(M_PI_2+(2-M_PI_2)*x));
  return multiplier*pow(k2Sqrt3*H(x)/(m_.q_*Rho(x)), m_.n_);
}

double Bound::Rho(double x) const{
  return ((ordinate_ == 'y') ? 1 : (1-x));
}

double Bound::dRho(double x) const{
  return ((ordinate_ == 'y') ? 0 : -1);
}

double Bound::S(double x) const{
  return ((ordinate_ == 'y') ? x : -1+Bound::kB+2*x); 
}

double Bound::dS(double x) const{
  return ((ordinate_ == 'y') ? 1 : 2); 
}

double Bound::Alpha(double x) const{
  return ((ordinate_ == 'y') ? M_PI_2 : M_PI_2/2); 
}

double Bound::dAlpha(double x) const{
  return 0;
}

double Bound::B1(double x) const{
  // std::cerr << x << ' '<< Rho(x)*dAlpha(x) << " alpha:" << Alpha(x)<< ' ' <<dRho(x) << ' ' << dS(x) << std::endl;
  return Rho(x)*dAlpha(x)+Alpha(x)*dRho(x)+dS(x);
}

double Bound::B2(double x) const{
  // std::cerr <<  x  << ' ' << Rho(x)*Alpha(x)+S(x) << std::endl; 
  
  // double q = 
  // std::cerr << x << ' ' << q << std::endl;
  // DCHECK(q>=0);
  return Rho(x)*Alpha(x)+S(x);;
}

double Bound::SigmaE(double x) const{
  // double h = H(x);
  // std::cerr << m_.q_ <<  ' ' << Rho(x) << ' ' << h << std::endl;
  return sqrt(3)/2*m_.q_*Rho(x)/H(x);
}

double Bound::H(double x) const{

  return (ordinate_ == 'y') ? m_.h1_/(2/M_PI*x+1) : m_.h1_/(1+(4/M_PI - 1)*x);//(1+(2/M_PI - 1)*x); //*exp(-Simpson::Integrate(0, x, kSimpsonStep, HFunctor(*this)));
}