#include "membrane.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <cassert>

#include "simpson.h"
#include "free_deformation.h"
#include "ideal_sliding.h"

using namespace std;

namespace {

bool IsNaN(double value){
  return value != value;
}

//TODO fix this const!
const double alphaElasticity = 0.226;


}  // namespace

Membrane::Membrane(double h0, double q, double n, double sigma_b, double a,double epsilon, int simpsonStep, int steps):
    h0_(h0), q_(q), n_(n), sigma_b_(sigma_b), a_(a), epsilon_(epsilon), simpsonStep_(simpsonStep), steps_(steps){
      da_ = (m_surface_.AlphaConstrained())/steps_;
      dx_ = (m_surface_.RightZero())/steps_;
      #ifdef __linux__
        num_threads_ = 1; //for debug constained step
        // num_threads_ = sysconf(_SC_NPROCESSORS_ONLN);
      #else
        num_threads_ = 2;
      #endif

      dstep_ = steps_/num_threads_;
}

double Membrane::operator()(double alpha) const{
  return (1/alpha-1/tan(alpha))*pow((2*h0_*sin(alpha)*sin(alpha)/(sqrt(3)*q_*alpha) -1), n_);
}


void Membrane::EasyIntegrate(int tid, vector<pair<double, double>> *v) {
  double t;
  int from = tid*dstep_;
  int to = (tid+1)*dstep_;

  FreeDeformation f(h0_, q_, n_);
  for(int i = from; i < to; i++){
    t = Simpson::Integrate(alphaElasticity+(i-1)*da_, alphaElasticity+(i+0)*da_, simpsonStep_, f);
    if(IsNaN(t))
      break;
    v->push_back(make_pair(t, i * da_));
  }
}

void Membrane::FreeStep(int steps){
  thread threads[num_threads_];
  vector<pair<double, double>> vectors[num_threads_];
  for(int i = 0; i<num_threads_; ++i)
    threads[i] = thread(bind(&Membrane::EasyIntegrate, this, i, &vectors[i]));

  for(int i = 0; i<num_threads_; ++i)
    threads[i].join();

  times_free_.clear();
  double offset = 0;
  for (int i = 0; i < num_threads_; ++i) {
    for (auto it = vectors[i].begin(); it != vectors[i].end(); ++it) {
      times_free_[it->first + offset] = it->second;
      offset += it->first;
    }
  }

  for(auto it = times_free_.cbegin(); it != times_free_.cend(); ++it) {
    cout<< it -> first << ' ' <<  it->second << ' ' << sin(it->second)/it->second*h0_ << endl;
  }
}

void Membrane::IntegrateConstrained(vector<pair<double, double>> *v){
  double t;
  int from = dstep_;
  int to = 1;
  IdealSliding is((*this));
  for(int i = from; i > to; i--){
    t = Simpson::Integrate((i-1)*dx_, (i+0)*dx_, 9, is);
    // cout << (i-1)*dx_ << ' ' << (i+0)*dx_ << endl;
    // assert(!IsNaN(t));
    // if(IsNaN(t))
    //   break;
    // cout << t << " " << i*dx_ << endl;
    // times_constrained_[t] = i*dx_;
    v->push_back(make_pair(t, i * dx_));
  }

}

void Membrane::ConstrainedStep(int steps){
  vector<pair<double, double>> vectors[num_threads_];
  
  // assert(1>2);
  Membrane::IntegrateConstrained(&vectors[0]);
  times_constrained_.clear();
    double offset = 0;
    for (int i = 0; i < num_threads_; ++i) {
      for (auto it = vectors[i].begin(); it != vectors[i].end(); ++it) {
        times_constrained_[it->first + offset] = it->second;
        offset += it->first;
      }
    }
}

void Membrane::IntegrateForAnimation(int steps){

  FreeStep(steps);

  // double dt_average, da_average;

  AverageDt(MeanValueDt());
  // da_average = MeanValueDa();

  // while(da_average > epsilon_){
  //   correctTimes();
  //   averageDt(meanValueDt());
  //   da_average = meanValueDa();
  // }
}

void Membrane::OutputResult(){
  // for(it = times_free_.begin(); it!=times_free_.end(); ++it){
  //   cout << it->first << ' ' << it->second << endl;
  // }
  // double x;
  IdealSliding is((*this));
  for(auto it = times_constrained_.cbegin(); it != times_constrained_.cend(); ++it) {
  // for(x=1; x>0; x-=0.1){
    // is.h(x);
    cout<< it -> first << ' ' << it->second << endl;
  }
}

double Membrane::ValueAsLine(double time, map<double, double>::iterator point1, map<double, double>::iterator point2){
  double k, b;
  k = (point2->second - point1->second)/(point2->first - point1->first);
  b = -k*point1->first+point1->second;

  return k*time + b;
}

void Membrane::CorrectTimes(){
  double alpha, tmp;
  map<double, double>::iterator it;

  // TODO (ukolshurika@): create  second map, fill it and swqap with times_free_.
  for(it = times_free_.begin(); it!=times_free_.end(); ++it){
    tmp = Simpson::Integrate(alphaElasticity, it->second, simpsonStep_, *this);
    alpha = it->second;
    // times_free_.erase(it);
    times_free_[tmp] = alpha;
  }
}

void Membrane::AverageDt(double dt_average){
  double delta = 0.0, next_time;

  map<double, double> new_times;
  map<double, double>::iterator next, cur;
  next = times_free_.begin();
  cur = next++;

  // TODO save DA & simpson values for correct data.
  while(next != times_free_.end()){
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


  times_free_.swap(new_times);
}

double Membrane::MeanValueDt(){
  assert(times_free_.size() > 2);

  double delta = 0.0;
  map<double, double>::iterator next, cur;
  next = times_free_.begin();
  cur = next++;

  while(next != times_free_.end()){
    delta += (next->first - cur->first);
    ++cur;
    ++next;
  }
  return delta/(times_free_.size()-1);
}

double Membrane::MeanValueDa(){
  return 0.0;
}
