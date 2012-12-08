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

struct LinearFunction {
  double operator () (double x) const {
    return x;
  }
};

struct SinusFunction {
  double operator () (double x) const {
    return sin(x);
  }
};

struct InverseSqrtFunction {
  double operator () (double x) const {
    return 1/sqrt(x);
  }
};

TEST(SimpsonTest, Constant) {
  ASSERT_TRUE(eq(0.0, Simpson::Integrate(0, 1, kNumSteps, ConstantFunction())));
}

TEST(SimpsonTest, Linear) {
  ASSERT_TRUE(eq(0.5, Simpson::Integrate(0, 1, kNumSteps, LinearFunction())));
}

TEST(SimpsonTest, Sinus) {
  ASSERT_TRUE(eq(2, Simpson::Integrate(0, M_PI, kNumSteps, SinusFunction())));
}

TEST(SimpsonTest, InverseSqrt) {
  ASSERT_TRUE(eq(2.0, Simpson::Integrate(1, 4, kNumSteps, InverseSqrtFunction())));
}