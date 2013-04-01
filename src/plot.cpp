#include <cmath>
#include <iostream>


#include "bound.h"
#include "membrane.h"
#include "simpson.h"

using namespace std;



int main(int argc, char **argv){
  double h0 = 2/100.0;
  double q = 2.65/88.3/1000;
  double n = 3.4;

  double t, x;

  Membrane m(q, h0, n);
  Bound b(m);
  if(argv[1][0] == 'd'){ // d means debug5
    for(x=1; x >= 0;x-=0.01){
      cout<< x << ' ' << (b.Rho(x)*b.dAlpha(x) + b.Alpha(x)*b.dRho(x) + b.dS(x))<< endl;
    }
  }else{
    while(cin >> t >> x){
      if(argv[1][0] == 'h')
        cout<< t/100000000.0 << ' ' << b.H(x)/h0 << endl;
      else if(argv[1][0] == 's')
        // b.SigmaE(x);
        cout<< t/100000000.0 << ' ' << b.SigmaE(x)*h0*1000+0.3 << endl;
      else if(argv[1][0] == 'a')
        cout<< t/100000000.0 << ' ' << sin(x)/(x) << endl;
      else if(argv[1][0] == 'e')
        cout<< t/100000000.0 << ' ' << q*x/h0/sin(x)/sin(x)*1000 << endl;


    }
  }
  return 0;
}
