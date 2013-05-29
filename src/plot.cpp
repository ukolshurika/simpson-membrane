#include <cmath>
#include <iostream>
#include <fstream>


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
  Bound b1(m);

  ifstream constrained_data_x("data/constrained_new.dat");
  ifstream free_data("data/free_new.dat");  

  if(argv[1][0] == 'd'){ // d means debug5
    for(x=1; x >= 0;x-=0.01){
      cout<< x << ' ' << b1.H(x) << endl;
    }
  }else{

    while(free_data >> t >> x){
       if(argv[1][0] == 'h')
        cout<< t/100000000.0 << ' ' << sin(x)/(x) << endl;
      else if(argv[1][0] == 's')
        cout<< t/100000000.0 << ' ' << q*x/sin(x)/sin(x)/h0*1000 << endl;
    }

    while(constrained_data_x >> t >> x){
      if(argv[1][0] == 'h')
        cout<< t/100000000.0 << ' ' << b1.H(x)/h0 << endl;
      else if(argv[1][0] == 's')
        cout << t/100000000.0 << ' ' << b1.SigmaE(x)*1000-92<< endl;
    }
  }

  return 0;
}