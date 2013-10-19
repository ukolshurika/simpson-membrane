#ifndef UTILS_H_
#define UTILS_H_

#include <cstdlib>
#include <iostream>

#if defined(DEBUG)
#define DCHECK(e) {             \
  if (!(e)) {               \
    std::cerr << "Check failed for " << (#e) << " "     \
              << "at " << __FILE__ << ":" << __LINE__ << std::endl; \
    *(int*) 0 = 0x31337;               \
  }                 \
}
#else
#define DCHECK(e)
#endif  // defined(DEBUG)

namespace utils{
  const double kEpsilon = 0.001;

  bool eql(double a, double b);
  bool IsNaN(double a);

  template<typename F>
  double KahanSum(int n, const F& f) {
    double s = 0, c = 0, t, y;
    for(int j=0; j<n; j++){
      y = f(j) - c;
      t = s + y;
      c = (t - s) -y;
      s = t;
    }
    return s;
  }
};

#endif  // UTILS_H_
