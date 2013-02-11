#include "membrane.h"

#include <iostream>
#include "ideal_sliding.h"

using namespace std;


int main(){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;
  double a = 1.0;
  double SigmaB = 88.3*1000000;
  double t, x;

  Membrane m(h0, q, n, SigmaB, a, 0.1, 999, 1);
  IdealSliding is(m.m_surface_, m);
  while(cin >> t >> x){
    cout<< t << ' ' << is.SigmaE(x) << endl;
  }
  return 0;
}
