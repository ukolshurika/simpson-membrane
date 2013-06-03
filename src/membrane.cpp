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
  cerr << h1_ << ' ' << h1_/h0_ << endl;
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

double p_k[10000], p_k1[10000];
double ds_k[10000], ds_k1[10000];
double delta_ds_k[10000], delta_ds_k1[10000];
double drho_k, drho_k1;
double sigma_k[10000], sigma_k1[10000];
double h_k[10000], h_k1[10000];
double dt = 1.5*100000, mu = 0.3;
ofstream constrained_data_h("data/constrained_h.dat");
ofstream constrained_data_s("data/constrained_s.dat");

void iteration(int iter){
  int i = 0;

  for(i = 1; i<=iter; ++i){
    delta_ds_k1[i] = pow((sigma_k1[i-1]+sigma_k[i])/(4/sqrt(3)-(sigma_k[i-1]+(sigma_k[i]))), n)*ds_k[i]*dt;
  }

  ds_k1[iter] = M_PI*(pow(k2Sqrt3*H(x)/(q_*Rho(x)), n_));

  for(i = 1; i<=iter; ++i){
    h_k1[i] = h_k[i](1-(pow(k2Sqrt3*H(x)/(q_*Rho(x)), m_.n_))));
  }

  sigma_k1[iter] = (m_.q_)/(m_.h0_*h_k1[iter]);

  for(i = iter-1; i>0; --i){
    sigma_k1[i] = sigma_k1[i+1]*h_k1[i+1]/h_k1[i] - mu*((ds_k1[i+1]*q_)/(h_k1[i])*h0_);
  }
}


void Membrane::constrained(int steps){
  int k;
  double x;

  for(i=10000; i>=0; i--){
    x=i/10000.0;
    sigma_k[i] = q_*(x)/sin(x)/sin(x)/h0_
  }

  for(k = 0; k< 10000; ++k){
    iteration(k);
    swap(delta_ds_k1, delta_ds_k);
    swap(h_k1, h_k);
    swap(sigma_k1, sigma_k);
    constrained_data_s << k*dt << " " << sigma_k[k] << endl;
    constrained_data_h << k*dt << " " << h_k[k] << endl;
  };

  for(k = 0; k< 10000; ++k){
    iteration(k);
    swap(delta_ds_k1, delta_ds_k);
    swap(h_k1, h_k);
    swap(sigma_k1, sigma_k);
    constrained_data_s << k*dt << " " << sigma_k[k] << endl;
    constrained_data_h << k*dt << " " << h_k[k] << endl;
  };
  // copy-paste part 4
}
