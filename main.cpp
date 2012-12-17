#include "membrane.h"

#include <iostream>
#include <ctime>

using namespace std;

int main(){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;

  int steps;
  cin >> steps;
  time_t start = time(NULL);
  Membrane m(h0, q, n, 0.1, 999, steps);
  m.IntegrateForAnimation(steps);
  time_t stop = time(NULL);

  cout << (stop - start) << endl;
  // m.OutputResult();
  return 0;
}