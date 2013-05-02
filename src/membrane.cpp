#include "membrane.h"

#include <cmath>
#include <iostream>
#include <sstream> 

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

void Membrane::constrained(int steps){
  
  double dx = Bound::kB/steps;
  double t;

  Bound b1(*this, 'y');
  Bound b2(*this, 'x');

  vector<pair<double, double>> v;
  for(double x = 0; x <= Bound::kB-1; x+=dx){
    t = Simpson::Integrate(x, x+dx, kSimpsonStep, b1);
    v.push_back(make_pair(t, x));
  }

  /*by y ordinate*/
  t_constrained_y_.clear();
  double offset = 0;
  double multiplire = sqrt(3)/2;
  double t_free_end = t_free_.back().first;
  double x_touch = Bound::kB - 1;
  double offset2 = 0;

  for (auto it = v.begin(); it != v.end(); ++it) {
    t_constrained_y_.push_back(make_pair(multiplire*(it->first + offset)+t_free_end, it->second));
    offset += it->first;
  }

  h1_ = b1.H(Bound::kB - 1);
  cerr << h1_ << endl;
  /*by x ordinate*/
  offset2 = Simpson::Integrate(0, x_touch, 999, b1);
  v.clear();
  for(double x = 0; x <= 1; x+=dx){
    t = Simpson::Integrate(x, x+dx, kSimpsonStep, b2);
    if(!utils::IsNaN(t))
    v.push_back(make_pair(t, x));
  }


  t_constrained_.clear();
  // offset2 = 0;
  
  for (auto it = v.begin(); it != v.end(); ++it) {
    t_constrained_.push_back(make_pair(multiplire*(it->first + offset2)+t_free_end, it->second));
    offset2 += it->first;
  }

 // for (auto it = t_constrained_.begin(); it != t_constrained_.end(); ++it) {
 //    t_constrained_y_.push_back(make_pair(it->first, it->second+Bound::kB));
 //    offset += it->first;
 //  }

}

/*double Membrane::ValueAsLine(double time, vector<pair<double, double>>::iterator point1, vector<pair<double, double>>::iterator point2){
  double k, b;
  k = (point2->second - point1->second)/(point2->first - point1->first);
  b = -k*point1->first+point1->second;
  return k*time + b;
}

double Membrane::MeanValueDt(){
  assert(t_free_.size() > 2);

  double delta = 0.0;
  vector<pair<double, double>>::iterator next, cur;
  next = t_free_.begin();
  cur = next++;

  while(next != t_free_.end()){
    delta += (next->first - cur->first);
    ++cur;
    ++next;
  }

  next = t_constrained_.begin();
  cur = next++;
  while(next != t_constrained_.end()){
    delta += (next->first - cur->first);
    ++cur;
    ++next;
  }

  next = t_constrained_y_.begin();
  cur = next++;
  while(next != t_constrained_y_.end()){
    delta += (next->first - cur->first);
    ++cur;
    ++next;
  }

  return delta/(t_free_.size() + t_constrained_y_.size() + t_constrained_.size() - 3);
}

void Membrane::averageDt(vector<pair<double, double>>* times, double dt_average){
  double delta, next_time;
  vector<pair<double, double>> new_times;
  vector<pair<double, double>><double, double>::iterator next, cur;
  next = (*times).begin();
  cur = next++;
  
  while(next != (*times).end()){
    delta = next->first - cur->first;

    if(new_times.size() == 0)
      new_times[cur->first] = new_times[cur->second];

    next_time = new_times.rbegin()->first + dt_average;

    if (delta >= dt_average){
      for(double i = new_times.rbegin()->first; i <= next->first; i += dt_average)
        new_times[i] = ValueAsLine(i, cur, next);

    }else if(next_time <= next->first){
      new_times[next_time] = ValueAsLine(next_time, cur, next);
    }

    ++cur;
    ++next;
  }


  (*times).swap(new_times);
}

void Membrane::averageAllDt(){
  double dt_average;
  dt_average = MeanValueDt();

  averageDt(t_free_,         dt_average);
  averageDt(t_constrained_,  dt_average);
  averageDt(t_constrained_y_, dt_average);
}*/

string Membrane::free_data_to_draw(double alpha){
  ostringstream data;
  data << "";
  #ifdef DRAWER
    data << " " << sin(alpha)/alpha*h0_ << " " << (1/sin(alpha)) << " " << (-1/tan(alpha));
  #endif
  return data.str();
}

string Membrane::constrained_y_data_to_draw(double x){
  ostringstream data;
  Bound b1(*this, 'y');
  data << "";
  #ifdef DRAWER
    data << " " << b1.H(x) << " " << b1.Rho(x) /*wide*/ << " " << x;
  #endif
  return data.str();
}

string Membrane::constrained_x_data_to_draw(double x){
  ostringstream data;
  Bound b2(*this, 'x');
  data << "";
  #ifdef DRAWER
    data << " " << b2.H(x) << " " << b2.Rho(x) << " " << 1-x;
  #endif
  return data.str();
}