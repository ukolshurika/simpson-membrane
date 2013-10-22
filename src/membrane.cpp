#include "membrane.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <cassert>

#include "bound.h"
#include "matrix.h"
#include "simpson.h"
#include "utils.h"

using namespace std;

namespace{

const double kSqrt3 = sqrt(3);
const int kSimpsonStep = 99;

double ValueAsLine(double time,
                   const pair<double, double>& point1,
                   const pair<double, double>& point2) {
  double k, b;
  k = (point2.second - point1.second)/(point2.first - point1.first);
  b = -k*point1.first+point1.second;
  return k*time + b;
}


void Kahan(double offset, const vector<pair<double, double>>& src, vector<pair<double, double>>* dst){
  DCHECK(dst->empty());
  for (auto& v : src)
    DCHECK(!utils::IsNaN(v.first));
  // cerr << 'v' <<endl;
  DCHECK(!utils::IsNaN(offset))
  double s = offset, c = 0, t, y;
  cerr << " #### " << s << endl;
  for(auto j=src.begin(); j!=src.end(); j++){
    y = j->first - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
    // DCHECK(kSqrt3*s/2.0 > offset);
    // (*dst).push_back(make_pair(kSqrt3*s/2.0, j->second));
    (*dst).push_back(make_pair(s, j->second));
  }
}

}  // namespace

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

Membrane::Membrane(double q, double h0, double n):q_(q), n_(n), h0_(h0) {
  alpha1_ = 0.41; // from Maxima flexible step(boolshit just get it from terraud)
  alpha2_ = 0.93;//atan(2*Bound::kB/(Bound::kB2));//M_PI/2;
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
  cerr << " ### "<< offset << endl;
  Kahan(offset, v, &t_free_);
  // for (auto it = v.begin(); it != v.end(); ++it) {
  //   t_free_.push_back(make_pair((it->first + offset), it->second));
  //   offset += it->first;
  // }
}

void Membrane::constrained(int steps){
  
  double dx = 1.0/steps;
  double t;

  Bound b1(*this, 'x');
  Bound b2(*this, 'y');

  vector<pair<double, double>> v;

  cerr<< dx << "QWEE   "<< Bound::kB << endl;
  for(double x = 0; x <= 1 - Bound::kB; x+=dx){
    t = Simpson::Integrate(x, x+dx, kSimpsonStep, b1); 
    DCHECK(t>0)
    if(!utils::IsNaN(t))
    v.push_back(make_pair(t, x));
  }

  //for(auto i : v) cerr << i.first << endl;

  /*by x ordinate*/
  t_constrained_.clear();
  
  double t_free_end = t_free_.back().first;
  cerr << "  !!!  " << t_free_end << endl;
  double offset2 = 0;

  DCHECK(!v.empty())

  Kahan(t_free_end, v, &t_constrained_);

  cerr << "  !!!  " << t_constrained_.front().first << endl;

  h1_ = b1.H(1-Bound::kB);
  cerr << h1_ << endl;
  /*by y ordinate*/
  
  offset2 = t_constrained_.back().first;//Simpson::Integrate(0, x_touch, 999, b1);
  
  v.clear();
  for(double x = 0; x+dx < 1; x+=dx){
    t = Simpson::Integrate(x, x+dx, kSimpsonStep, b2);
    DCHECK(t>0)
    if(!utils::IsNaN(t))
      v.push_back(make_pair(t, x));
  }

  t_constrained_y_.clear();
  cerr<< "  @@@  "<< offset2 <<endl;
  Kahan(offset2, v, &t_constrained_y_);
}

double Membrane::MeanValueDt(){
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
  assert((*times).size() > 2);
  vector<pair<double, double>> new_times;
  vector<pair<double, double>>::iterator next, cur;
  next = (*times).begin();
  cur = next++;

  while(next != (*times).end()){
    delta = next->first - cur->first;

    if(new_times.size() == 0)
      new_times.push_back(make_pair(cur->first,cur->second));

    next_time = new_times.rbegin()->first + dt_average;

    if (delta >= dt_average){
      for(double i = new_times.rbegin()->first; i <= next->first; i += dt_average)
        new_times.push_back(make_pair(i, ValueAsLine(i, *cur, *next)));

    }else if(next_time <= next->first){
      new_times.push_back(make_pair(next_time,
                                    ValueAsLine(next_time, *cur, *next)));
    }
    ++cur;
    ++next;
  }

  (*times).swap(new_times);
}

void Membrane::averageAllDt(){
  double dt_average;
  dt_average = MeanValueDt();

  averageDt(&t_free_,         dt_average);
  averageDt(&t_constrained_,  dt_average);
  averageDt(&t_constrained_y_, dt_average);
}

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