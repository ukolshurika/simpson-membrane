#include "membrane.h"

#include <iostream>

using namespace std;

int main(){
  double h0 = 0.02;
  double q = 2650.0/88.3/1000000;
  double n = 3.4;

  int steps;
  cin >> steps;
  Membrane m(h0, q, n, 0.1, 999, steps);
  m.IntegrateForAnimation(steps);
  m.OutputResult();
  return 0;
}