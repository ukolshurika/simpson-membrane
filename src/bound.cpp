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
    if(b.ordinate_=='x'){
      DCHECK(x>=0);
      return b.B1(x)/b.B2(x);
    }else{
      return (2-M_PI_2)/(Bound::kB+(M_PI_2-1))+(2-M_PI_2)*x;
    }
  }

  Bound b;
};

const double Bound::kB = 0.5;
const double Bound::kB2 = 0.25;

Bound::Bound(const Membrane& m, char ordinate):m_(m), ordinate_(ordinate){};

double Bound::operator()(double x) const{
  // std::cerr << "CONSTRAINED" << std::endl;
   // std::cerr << "!!!!"<< B1(x) << ' ' << B2(x) << std::endl;
  if(utils::IsNaN(B1(x)/B2(x))) std::cerr << x << std::endl;
  DCHECK(B1(x)>0);
  DCHECK(B2(x)>0);
  DCHECK(x>=0);
  DCHECK(x<=1);
  double multiplier;
  multiplier = (ordinate_ == 'x') ? B1(x)/B2(x) : ((2-M_PI_2)/(M_PI_2+(2-M_PI_2)*x));
  return multiplier*pow(k2Sqrt3*H(x)/(m_.q_*Rho(x)), m_.n_);
}

double Bound::Rho(double x) const{
  return ((ordinate_ == 'x') ? (kB*kB+(1-x)*(1-x))/2/kB : (1-x));
}

double Bound::dRho(double x) const{
  return ((ordinate_ == 'x') ? (2*(x-1)) : -1);
}

double Bound::S(double x) const{
  return ((ordinate_ == 'x') ? x : 1-Bound::kB+2*x); 
}

double Bound::dS(double x) const{
  return ((ordinate_ == 'x') ? 1 : 2); 
}

double Bound::Alpha(double x) const{
  // std::cerr << "alpha: " << (2*kB*(1-x))/((1-x)*(1-x)+kB2) << " x:" << x << std::endl;
  return ((ordinate_ == 'x') ? atan((2*kB*(1-x))/((1-x)*(1-x)+kB2)) : M_PI_2/2); 
}

double Bound::dAlpha(double x) const{
  double x2 = x*x;
  double x3 = x*x*x;
  double x4 = x3*x;
  return (2*kB*(x-kB-1)*(x+kB-1))/(x4 - 4*x3 + 6*kB2*x2 + 6*x2 - 12*kB2*x - 4*x + kB2*kB2 + 6*kB2 + 1);
}

double Bound::B1(double x) const{
  DCHECK(dS(x)>0);
  return Alpha(x)*dRho(x)+dAlpha(x)*Rho(x)+dS(x);
}

double Bound::B2(double x) const{
  DCHECK(x<=1);
  DCHECK(Rho(x)>0);
  DCHECK(Alpha(x)>0);
  DCHECK(S(x)>=0);
  return Rho(x)*Alpha(x)+S(x);;
}

double Bound::SigmaE(double x) const{
  // double h = H(x);
  // std::cerr << m_.q_ <<  ' ' << Rho(x) << ' ' << h << std::endl;
  return sqrt(3)/2*m_.q_*Rho(x)/H(x);
}

double Bound::H(double x) const{
  DCHECK(x>=0);
  return (ordinate_ == 'x') ? m_.h1_*exp(Simpson::Integrate(x, 1-kB, kSimpsonStep, HFunctor(*this))) : m_.h1_/(1+(4/M_PI - 1)*x);//(1+(2/M_PI - 1)*x); //*exp(-Simpson::Integrate(0, x, kSimpsonStep, HFunctor(*this)));
} 