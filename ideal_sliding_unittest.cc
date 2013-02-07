#include "ideal_sliding.h"

#include <iostream>
#include <cmath>

#include "gtest/gtest.h"
#include "membrane.h"
#include "matrix_surface.h"
#include "free_deformation.h"

namespace {
const double kEpsilon = 1e-9;

bool IsNaN(double value){
  return value != value;
}

bool eql(double u, double v){
  return abs(u - v) < kEpsilon;
}

Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
MatrixSurface ms;
IdealSliding is(ms, m, 0.2);

}


TEST(Circle, center00){
  ASSERT_TRUE(eql(0, Circle::center(1, M_PI)));
}

TEST(Circle, 2centerm10){
  ASSERT_TRUE(eql(-1, Circle::center(1, 0)));
}

TEST(Circle, centerm10){
  ASSERT_TRUE(eql(-1, Circle::center(2, M_PI/2)));
}

TEST(IdealSlidingTest, Alpha) {
  double x;
  for(x = 0; x<1; x+=0.1){
    ASSERT_TRUE(!IsNaN(is.Alpha(x)));
    ASSERT_TRUE(is.Alpha(x) >= 0);
  }
}

TEST(IdealSlidingTest, dAlpha) {
  double x;
  for(x = 0; x<1; x+=0.1)
    ASSERT_TRUE(!IsNaN(is.dAlpha(x)));
}

TEST(IdealSlidingTest, Rho) {
  double x;
  for(x = 0; x<0.99; x+=0.1){
    ASSERT_TRUE(!IsNaN(is.Rho(x)));
    ASSERT_TRUE(is.Rho(x) >= 0);
  }
}

TEST(IdealSlidingTest, dRho) {
  double x;
  for(x = 0; x<1; x+=0.1)
    ASSERT_TRUE(!IsNaN(is.dRho(x)));
}

TEST(IdealSlidingTest, S){
  double x;
  for(x = 0; x<1; x+=0.1){
    ASSERT_TRUE(!IsNaN(is.S(x)));
    ASSERT_TRUE(is.S(x) >= 0);
  }
}

TEST(IdealSlidingTest, dS){
  double x;
  for(x = 1; x>0; x-=0.1){
    ASSERT_TRUE(!IsNaN(is.dS(x)));
    ASSERT_TRUE(is.dS(x) >= 0);
  }
}

TEST(IdealSlidingTest, Operator){
  double x;
  for(x = 0.01; x<1; x+=0.1){
    ASSERT_TRUE(!IsNaN(is(x)));
  }
}

TEST(IdealSlidingTest, B1){
  double x;
  for(x = 0; x<1; x+=0.1)
    ASSERT_TRUE(!IsNaN(is.B1(x)));
}

TEST(IdealSlidingTest, B2){
  double x;
  for(x = 0; x<1; x+=0.1)
    ASSERT_TRUE(!IsNaN(is.B2(x)));
}

TEST(IdealSlidingTest, h){
  double x;
  for(x = 0; x<1; x+=0.1){
    ASSERT_TRUE(!IsNaN(is.h(x)));
    ASSERT_TRUE(is.h(x)>0);
  }
}

