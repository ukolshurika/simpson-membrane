#include "matrix_surface.h"

#include <iostream>
#include <cmath>

#include "gtest/gtest.h"



namespace {
const double kEpsilon = 1e-9;

bool IsNaN(double value){
  return value != value;
}

bool eql(double u, double v){
  return abs(u - v) < kEpsilon;
}

MatrixSurface ms;
}

TEST(MatrixSurfaceTest, Operator) {
  double x = 0;
  for(x = 0; x< 1; x+=0.1)
    ASSERT_TRUE(eql(ms(x), -x*x+1));
}

TEST(MatrixSurfaceTest, AlphaConstrained) {
  ASSERT_TRUE(eql(M_PI-atan(-2), ms.AlphaConstrained()));
}

TEST(MatrixSurfaceTest, SecondDerivative) {
  double x = 0;
  for(x = 0; x < 1; x+=0.1)
    ASSERT_TRUE(eql(-2, ms.SecondDerivative(x)));
}

TEST(MatrixSurfaceTest, Derivative) {
  double x = 0;
  for(x = 0; x < 1; x+=0.1)
    ASSERT_TRUE(eql(-2*x, ms.Derivative(x)));
}


TEST(MatrixSurfaceTest, dNormal) {
  double x;
  for(x=0.1; x < 1; x+=0.1)
    ASSERT_TRUE(eql(-1/2.0/x, ms.dNormal(x)));

}


TEST(MatrixSurfaceTest, RightZero) {
  ASSERT_TRUE(eql(1, ms.RightZero()));
}