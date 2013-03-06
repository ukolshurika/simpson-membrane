#include "membrane.h"

#include <iostream>
#include "ideal_sliding.h"
#include "simpson.h"
using namespace std;



int main(int argc, char **argv){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;
  double a = 1.0;
  double SigmaB = 88.3*1000000;
  double t, x;

  Membrane m(h0, q, n, SigmaB, a, 0.1, 999, 1);
  IdealSliding is(m);
  if(argv[1][0] == 'd'){ // d means debug
    for(x=1; x >= 0;x-=0.01){
      t = is.h(x);
      cout<< x << ' ' << is.Rho(x) << endl;
    }
  }else{
    while(cin >> t >> x){
      if(argv[1][0] == 'h')
        cout<< t << ' ' << is.h(x) << endl;
      else if(argv[1][0] == 's')
        cout<< t << ' ' << is.SigmaE(x) << endl;
    }
  }
  return 0;
}
