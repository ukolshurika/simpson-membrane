#include "membrane.h"

#include <cmath>
#include <iostream>

#include "bound.h"
#include "matrix.h"
#include "simpson.h"
#include "utils.h"

using namespace std;

namespace{
  const double kSqrt3 = sqrt(3);
  const int kSimpsonStep = 99;
}

struct Free {
  Free(double h0, double q, double n):q_(q), h0_(h0), n_(n){};
  double operator()(double alpha) const{
    // cerr << "FREE" << endl;
    // cerr<< (1/alpha-1/tan(alpha)) << endl;
    return (1/alpha-1/tan(alpha))*pow(((2*h0_*sin(alpha)*sin(alpha))/(kSqrt3*q_*alpha)-1), n_);
  }

  double h(double alpha){
    return sin(alpha)/alpha*h0_;
  }

  private:
  double q_;
  double h0_;
  double n_;
};

Membrane::Membrane(double q, double h0, double n):q_(q), h0_(h0), n_(n){
  alpha1_ = 0.41; // from Maxima flexible step(boolsh it just get it from terraud)
  alpha2_ = M_PI/2;
  h1_ = sin(alpha2_)/alpha2_*h0_;
  cerr << h1_ << endl;
}

void Membrane::free(int steps){
  double dalpha = (alpha2_ - alpha1_)/steps;
  double t;
  Free f(h0_, q_, n_);

  vector<pair<double, double>> v;

  for(double a = alpha1_; a < alpha2_; a+=dalpha){
    t = Simpson::Integrate(a, a+dalpha, kSimpsonStep, f);
    v.push_back(make_pair(t, a));
  }

  t_free_.clear();
  double offset = 0;

  for (auto it = v.begin(); it != v.end(); ++it) {
    t_free_.push_back(make_pair((it->first + offset), it->second));
    offset += it->first;
  }
}

void Membrane::constrained(int steps){
  
  double dx = Bound::kB/steps;
  double t;

  Bound b1(*this, 'y');
  Bound b2(*this, 'x');

  vector<pair<double, double>> v;
  for(double x = 0; x <= Bound::kB; x+=dx){
    t = Simpson::Integrate(x, x+dx, kSimpsonStep, b1);
    v.push_back(make_pair(t, x));
  }

  /*by y ordinate*/
  t_constrained_y_.clear();
  double offset = 0;
  double multiplire = sqrt(3)/2;
  double t_free_end = t_free_.back().first;

  for (auto it = v.begin(); it != v.end(); ++it) {
    t_constrained_y_.push_back(make_pair(multiplire*(it->first + offset)+t_free_end, it->second));
    offset += it->first;
    // cerr << it->first << ' ' << it->second << endl;
  }

  h1_ = b1.H(Bound::kB - 1);

  /*by x ordinate*/
  v.clear();
  for(double x = 0; x <= 1; x+=dx){
    t = Simpson::Integrate(x, x+dx, kSimpsonStep, b2);
    if(!utils::IsNaN(t))
    v.push_back(make_pair(t, x));
  }


  t_constrained_.clear();
  // offset = 0;
  
  for (auto it = v.begin(); it != v.end(); ++it) {
    t_constrained_.push_back(make_pair(multiplire*(it->first + offset)+t_free_end, it->second));
    offset += it->first;
  }



}