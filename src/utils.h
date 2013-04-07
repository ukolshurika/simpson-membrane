#ifndef UTILS_H_
#define UTILS_H_

#include <cstdlib>
#include <iostream>

#if defined(DEBUG)
#define DCHECK(e) {             \
  if (!(e)) {               \
    std::cerr << "Check failed for " << (#e) << " "     \
              << "at " << __FILE__ << ":" << __LINE__ << std::endl; \
    exit(-1);               \
  }                 \
}
#else
#define DCHECK(e)
#endif  // defined(DEBUG)

namespace util{
  const double kEpsilon = 0.001;

  bool eql(double a, double b);
  bool IsNan(double a);

};

#endif  // UTILS_H_
