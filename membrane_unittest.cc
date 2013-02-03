#include "membrane.h"

#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

const double kEpsilon = 1e-9;
const int kNumSteps = 999;


bool eql(double u, double v) {
  return abs(u - v) < kEpsilon;
}

//===============ValueASLine==================//
TEST(MembraneTest, ConstantValueAsLine) {
  std::map<double, double> points;
  points[0] = 0; points[1] = 0;
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  ASSERT_TRUE(eql(0.0, m.ValueAsLine(0.5, points.begin(), points.end())));
}

TEST(MembraneTest, PositiveKValueAsLine) {
  std::map<double, double> points;
  points[0] = 0; points[1] = 1;
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  ASSERT_TRUE(eql(0.5, m.ValueAsLine(0.5, points.begin(), points.end())));
}

TEST(MembraneTest, NegativeKValueAsLine) {
  std::map<double, double> points;
  points[0] = 0; points[1] = -1;
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  ASSERT_TRUE(eql(-0.5, m.ValueAsLine(0.5, points.begin(), points.end())));
}

TEST(MembraneTest, positiveKValueAsLine) {
  std::map<double, double> points;
  points[0] = 2; points[3] = 7;
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  ASSERT_TRUE(eql(1.0, m.ValueAsLine(0.5, points.begin(), points.end())));
}

//===============MeanValueDt==================//
TEST(MembraneTest, EqualMeanValueDt) {
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  for(double i = 0; i<10; ++i)
    m.times_free_[i] = 1.0;
  ASSERT_TRUE(eql(1.0, m.MeanValueDt()));
}

TEST(MembraneTest, LinearMeanValueDt) {
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  for(double i = 0; i<10; ++i)
    m.times_free_[i*2] = 1.0;

  ASSERT_TRUE(eql(2.0, m.MeanValueDt()));
}

TEST(MembraneTest, QuadrickMeanValueDt) {
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  for(double i = 0; i<10; ++i)
    m.times_free_[i*i] = 1.0;

  ASSERT_TRUE(eql(9.0, m.MeanValueDt()));
}


//===============AverageDt==================//
TEST(MembraneTest, EqualAverageDt) {
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  for(double i = 0; i<10; ++i)
    m.times_free_[i] = i*2;

  m.AverageDt(m.MeanValueDt());

  for(double i = 0; i<10; ++i)
    ASSERT_TRUE(eql(i*2, m.times_free_[i]));
}

TEST(MembraneTest, LinearAverageDt) {
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  for(double i = 0; i<10; ++i)
    m.times_free_[i*2] = i*0.5+4;

  m.AverageDt(m.MeanValueDt());

  for(double i = 0; i<10; ++i)
    ASSERT_TRUE(eql(i*0.5+4, m.times_free_[i*2]));
}

TEST(MembraneTest, QuadAverageDt) {
  Membrane m(0.02, 2650.0/88.3/1000000, 3.4, 88.3*1000000, 1.0,  0.1, 999, 1000);
  for(double i = 0; i<10; ++i)
    m.times_free_[i*i] = i;

  m.AverageDt(m.MeanValueDt());

  // std::map<double, double>::iterator i;
  // for(i=m.times_free_.begin(); i!=m.times_free_.end(); ++i)
  //   std::cout << i->first << std::endl;

  for(double i = 0; i<81; i+=9)
    ASSERT_TRUE(eql(sqrt(i), m.times_free_[i]));
}