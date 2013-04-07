#ifndef MATRIX_H_
#define MATRIX_H_

class Matrix{
public:
  double operator()(double x);
  static double Derivative(double x);
  static double RZero();

  const static double kK;
  const static double kB;
  const static double kKBSquare;
};


#endif //MATRIX_H_