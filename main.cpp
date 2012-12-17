#include "membrane.h"

#include <iostream>
#include <ctime>

using namespace std;

int main(){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;
  
  std::clock_t start;
  int steps;
  cin >> steps;


  start = clock();
  
  Membrane m(h0, q, n, 0.1, 999, steps);
  m.IntegrateForAnimation(steps);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);

  cout << (( std::clock() - start ) / (double) CLOCKS_PER_SEC) << endl;
  // m.OutputResult();
  return 0;
}