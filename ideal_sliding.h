#ifndef IDEAL_SLIDING_H_
#define IDEAL_SLIDING_H_

#include <map>
#include <vector>
#include <utility>
#include <cmath>

#include "matrix_surface.h"
#include "membrane.h"

class IdealSliding{
  public:
  IdealSliding(const MatrixSurface& ms, const Membrane& m);

  double operator () (double x)const;
  double Alpha(double x) const;
  double dAlpha(double x) const;
  double S(double x) const;
  double dS(double x) const;
  double Rho(double x) const;
  double dRho(double x) const;
  double B1(double x) const;
  double B2(double x) const;
  double h(double x) const;
  double SigmaE(double x) const;

  MatrixSurface ms_;
  Membrane m_;
  double h1_;
};

struct Circle{
  static double center(double a, double alpha);
};

struct HFunctor{
  HFunctor(const IdealSliding& is); 
  double operator () (double x) const;
  IdealSliding is_;
};

struct SFunctor{
  SFunctor(const IdealSliding& is); 
  double operator () (double x) const;
  IdealSliding is_;
};

#endif  // IDEAL_SLIDING_H