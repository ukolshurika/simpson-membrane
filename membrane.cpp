#include "membrane.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <map>

#include "simpson.h"

using namespace std;

namespace {

bool IsNaN(double value){
  return value != value;
}

//TODO fix this const!
const double alphaStart = 0.226;

}  // namespace

Membrane::Membrane(double h0, double q, double n, double epsilon, int simpsonStep, int steps):
    h0_(h0), q_(q), n_(n), epsilon_(epsilon), simpsonStep_(simpsonStep), steps_(steps){
      da_ = (M_PI)/steps_;
      num_threads_ = 2;
      dstep_ = steps_/num_threads_;
}

double Membrane::operator()(double alpha) const{
  return (1/alpha-1/tan(alpha))*pow((2*h0_*sin(alpha)*sin(alpha)/(sqrt(3)*q_*alpha) -1), n_);
}

void Membrane::EasyIntegrate(int tid, map<double, double>* m){
  double t = 0.0, tmp;
  int from = tid*dstep_;
  int to = (tid+1)*dstep_;

  for(int i = from; i < to; i++){
    tmp = Simpson::Integrate(alphaStart+(i-1)*da_, alphaStart+(i+0)*da_, simpsonStep_, *this);
    if(IsNaN(tmp))
      break;
    // t += tmp;
    (*m)[tmp] = i*da_;
  }
}

void Membrane::IntegrateForAnimation(int steps){

  thread threads[num_threads_];
  map<double, double> maps[num_threads_];
  for(int i = 0; i<num_threads_; ++i)
    threads[i] = thread(bind(&Membrane::EasyIntegrate, this, i, &maps[i]));

  for(int i = 0; i<num_threads_; ++i)
    threads[i].join();

  double offset = 0, last = 0;
  for (int i = 0; i < num_threads_; ++i) {
    for (auto it = maps[i].begin(); it != maps[i].end(); ++it) {
      times_[it->first + offset] = it->second;
      offset+=it->first;
      // last = it->first;
    }
    // offset = last;
  }
  // EasyIntegrate(steps); // get initial values

  double dt_average, da_average;
  // dt_average = 0.0;
  // da_average = 1.0;


  AverageDt(MeanValueDt());
  da_average = MeanValueDa();
  cout << MeanValueDa();
  // while(da_average > epsilon_){
  //   correctTimes();
  //   averageDt(meanValueDt());
  //   da_average = meanValueDa();
  // }
}

void Membrane::OutputResult(){
  map<double, double>::iterator it;

  for(it = times_.begin(); it!=times_.end(); ++it){
    cout << it->first << ' ' << it->second << endl;
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

  // TODO (ukolshurika@): create  second map, fill it and swqap with times_.
  for(it = times_.begin(); it!=times_.end(); ++it){
    tmp = Simpson::Integrate(alphaStart, it->second, simpsonStep_, *this);
    alpha = it->second;
    // times_.erase(it);
    times_[tmp] = alpha;
  }
}

void Membrane::AverageDt(double dt_average){
  double delta = 0.0, next_time;

  map<double, double> new_times;
  map<double, double>::iterator next, cur;
  next = times_.begin();
  cur = next++;

  // TODO save DA & simpson values for correct data.
  while(next != times_.end()){
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

   
  times_.swap(new_times);
}

double Membrane::MeanValueDt(){
  assert(times_.size() > 2);

  double delta = 0.0;
  map<double, double>::iterator next, cur;
  next = times_.begin();
  cur = next++;

  while(next != times_.end()){
    delta += (next->first - cur->first);
    ++cur;
    ++next;
  }
  return delta/(times_.size()-1);
}

double Membrane::MeanValueDa(){

  return 0.0;
}
