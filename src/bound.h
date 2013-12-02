#ifndef BOUND_H_
#define BOUND_H_

#include "membrane.h"


class Bound{
public:
  const static int kSimpsonStep = 99;
  const static double kB;
  const static double kB2;

  Bound(const Membrane& m, char ordinate);

  double operator()(double x) const;
  double B1(double x) const;
  double B2(double x) const;

  double Rho(double x) const;
  double dRho(double x) const;
  double S(double x) const;
  double dS(double x) const;
  double Alpha(double x) const;
  double dAlpha(double x) const;
  double q(double l) const;

  double SigmaE(double x) const;
  double H(double x) const;

  Membrane m_;
  char ordinate_;
};

#endif //BOUND_H_