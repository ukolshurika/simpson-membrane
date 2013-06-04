#include "membrane.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>

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

#define N 2
const double k2Sqrt3 = 2/sqrt(3);
double p_k[N], p_k1[N];
double ds_k[N], ds_k1[N];
double delta_ds_k[N], delta_ds_k1[N];
double drho_k[N], drho_k1[N];
double rho_k[N], rho_k1[N];
double sigma_k[N], sigma_k1[N];
double h_k[N], h_k1[N];
double dt = 1.5*100000, mu = 0.3;
ofstream constrained_data_h("data/constrained_h.dat");
ofstream constrained_data_s("data/constrained_s.dat");

void Membrane::iteration_vert(int iter){
  int i = 0;

  for(i = 1; i<N; ++i){
    delta_ds_k1[i] = pow((sigma_k[i-1]+sigma_k[i])/(4/sqrt(3)-(sigma_k[i-1]+(sigma_k[i]))), n_)*ds_k[i]*dt;
  }

  cout << endl;

  ds_k1[iter] = M_PI*(pow(k2Sqrt3*h_k[i]/(q_), n_));

  for(i = 1; i<=iter; ++i){
    h_k1[i] = h_k[i]*(1-(pow(k2Sqrt3*h_k[i]/(q_), n_)));
  }

  sigma_k1[iter] = (q_)/(h0_*h_k1[iter]);

  for(i = iter-1; i>0; --i){
    sigma_k1[i] = sigma_k1[i+1]*h_k1[i+1]/h_k1[i] - mu*((ds_k1[i+1]*q_)/(h_k1[i])*h0_);
  }
}

void Membrane::iteration_gorizontal(int iter){
  int i = 0;
  double sum = 0;

  for(i = 1; i<=iter; ++i){
    delta_ds_k1[i] = pow((sigma_k[i-1]+sigma_k[i])/(4/sqrt(3)-(sigma_k[i-1]+(sigma_k[i]))), n_)*ds_k[i]*dt;
    sum += delta_ds_k1[i];
    cout << sum << " ";
  }
  cout << endl;
  ds_k1[iter] = (sum+M_PI_2*(pow(k2Sqrt3*h_k[iter]/(q_*rho_k[iter]), n_))*rho_k[iter])/4;

  for(i = 0; i<=iter; ++i){
    drho_k1[i] = -(sum + 2*ds_k1[iter]);
  }

  sum = 0;
  for(i = 0; i<=iter; ++i){
    sum+=ds_k1[i];
    rho_k1[i] = 1-(sum);
  }

  cout << "H: ";
  for(i = 0; i<=iter; ++i){
    h_k1[i] = h_k[i]*(1-(pow(k2Sqrt3*h_k[i]/(q_*rho_k[i]), n_)));
    cout << h_k[i]/(q_*rho_k[i])<< " ";
  }



  sigma_k1[iter] = (q_)/(h0_*h_k1[iter]);
  cout <<  h_k1[iter] << endl;

  for(i = iter-1; i>=0; --i){
    sigma_k1[i] = sigma_k1[i+1]*h_k1[i+1]/h_k1[i] - mu*((ds_k1[i+1]*q_)/(h_k1[i])*h0_);
    cout << i<< ":" << sigma_k1[i+1] << " ";
  }
  cout << endl;
}

void Membrane::constrained(int steps){
  int k, i;
  double x, tmp;

  memset(p_k, 0, N*sizeof(double));        memset(p_k1, 0, N*sizeof(double));
  memset(rho_k, 0, N*sizeof(double));      memset(rho_k1, 0, N*sizeof(double));
  memset(drho_k, 0, N*sizeof(double));     memset(drho_k1, 0, N*sizeof(double));
  memset(ds_k, 0, N*sizeof(double));       memset(ds_k1, 0, N*sizeof(double));
  memset(delta_ds_k, 0, N*sizeof(double)); memset(delta_ds_k1, 0, N*sizeof(double));
  memset(sigma_k, 0, N*sizeof(double));    memset(sigma_k1, 0, N*sizeof(double));
  memset(h_k, 0, N*sizeof(double));        memset(h_k1, 0, N*sizeof(double));

  for(i=0; i<N; ++i){
    sigma_k[i] = q_/sin(1)/sin(1)/h0_;
    h_k[i] = h1_;
    cout << h_k[i] << " !! ";
  }

  cout << endl << " jjjjj " << endl;

  for(k = 2; k<N; ++k){
    iteration_vert(k);
    cout << k<< endl << "================"<< endl;
    for(i=0; i<N; ++i){
     // cout << i<< ":" << sigma_k[i] << " ";
    }
    cout << endl;

    
    swap(delta_ds_k1, delta_ds_k);
    swap(h_k1, h_k);
    swap(sigma_k1, sigma_k);
    
    constrained_data_s << k*dt << " " << sigma_k[k] << endl;
    constrained_data_h << k*dt << " " << h_k[k] << endl;
  };

  tmp = sigma_k[k];
  for(i=0; i<N; ++i){
    sigma_k[i] = tmp;
  }

  for(k = 0; k< N; ++k){
    iteration_gorizontal(k);
    swap(delta_ds_k1, delta_ds_k);
    swap(h_k1, h_k);
    swap(sigma_k1, sigma_k);
    swap(drho_k1, drho_k);
    swap(rho_k1, rho_k);
    constrained_data_s << k*dt + N*dt<< " " << sigma_k[k] << endl;
    constrained_data_h << k*dt + N*dt<< " " << h_k[k] << endl;
  };

}
