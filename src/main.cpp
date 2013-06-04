#include <iostream>
#include <fstream>
#include <cmath>

#include "membrane.h"
#include "bound.h"

using namespace std;


double h(double alpha){
  return sin(alpha)/alpha*0.02;
}

int main(){
  double h0 = 2/100.0;
  double q = 2.65*1000/88.3/1000000;
  double n = 3.4;

  Membrane m(q, h0, n);
  ofstream free_data("data/free_new.dat");
  ofstream constrained_data_x("data/constrained_x.dat");
  ofstream constrained_data_y("data/constrained_y.dat");

  m.free(999);
  m.constrained(999);


  for(auto i = m.t_free_.cbegin(); i != m.t_free_.cend(); ++i){
    free_data<< i -> first << ' '<< i->second << endl;
  }

  // for(auto i = m.t_constrained_.cbegin(); i != m.t_constrained_.cend(); ++i){
  //   constrained_data_x<< i -> first << ' ' << i->second << endl;
  // }


  // for(auto i = m.t_constrained_y_.cbegin(); i != m.t_constrained_y_.cend(); ++i){
  //   constrained_data_y<< i -> first << ' ' << i->second << endl;
  // }

  cerr << "C0MPLETE" << endl;
  return 0;
}