#include "membrane.h"

#include <iostream>
#include <ctime>
#include <sys/time.h>

using namespace std;

//TODO get correct value from article


int main(){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;
  double a = 1.0;
  double SigmaB = 1.0;//88.3*1000000;

  // timeval tim;
  int steps;
  cin >> steps;

  // gettimeofday(&tim,NULL);
  // double t1 = tim.tv_sec+(tim.tv_usec/1000000.0);

  Membrane m(h0, q, n, SigmaB, a, 0.1, 999, steps);
  m.ConstrainedStep(steps);
  // m.IntegrateForAnimation(steps);

  // gettimeofday(&tim,NULL);
  // double t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
  // cout << t2-t1 << endl;
  m.OutputResult();
  return 0;
}
