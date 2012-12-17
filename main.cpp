#include "membrane.h"

#include <iostream>
#include <ctime>

using namespace std;

int main(){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;
  
  timespec start, stop;
  int steps;
  cin >> steps;


  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  
  Membrane m(h0, q, n, 0.1, 999, steps);
  m.IntegrateForAnimation(steps);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);

  cout << ((stop.tv_nsec)/100000000.0 - (start.tv_nsec)/100000000.0) << endl;
  // m.OutputResult();
  return 0;
}