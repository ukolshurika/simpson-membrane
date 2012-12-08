#include "simpson.h"

#include <cmath>

#include "gtest/gtest.h"

const double kEpsilon = 1e-9;
const int kNumSteps = 999;

bool eq(double u, double v) {
  return abs(u - v) < kEpsilon;
}

// Constant function that always returns 0.
struct ConstantFunction {
  double operator () (double) const {
    return 0;
  }
};

TEST(SimpsonTest, Constant) {
  ASSERT_TRUE(eq(0.0, Simpson::Integrate(0, 1, kNumSteps, ConstantFunction())));
}
