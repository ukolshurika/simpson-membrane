#ifndef BOUND_H_
#define BOUND_H_

#include "membrane.h"


class Bound{
public:


const static int kSimpsonStep = 99;

  Bound(const Membrane& m);

  double operator()(double x) const;
  double B1(double x) const;
  double B2(double x) const;

  double Rho(double x) const;
  double dRho(double x) const;
  double S(double x) const;
  double dS(double x) const;
  double Alpha(double x) const;
  double dAlpha(double x) const;


  double SigmaE(double x) const;
  double H(double x) const;

  void PrintX0X1(int i);

  Membrane m_;
};

#endif //BOUND_H_