#ifndef UTILS_H_
#define UTILS_H_

#include <cstdlib>
#include <iostream>

#if defined(ALL_CHECKS)
#define CHECK(e) {							\
  if (!(e)) {								\
    std::cerr << "Check failed for " << (#e) << " "			\
              << "at " << __FILE__ << ":" << __LINE__ << std::endl;	\
    exit(-1);								\
  }									\
}
#else
#define CHECK(e)
#endif  // defined(ALL_CHECKS)

#endif  // UTILS_H_
