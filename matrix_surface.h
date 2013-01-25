#ifndef MATRIX_SURFACE_H_
#define MATRIX_SURFACE_H_

#define DELTA 0.001
class MatrixSurface{
  public:
    double operator () (double x) const;
    double Derivative (double x) const;
    double SecondDerivative (double x) const;
    double AlphaConstrained () const;
    double RightZero () const;
    double Normal(double x) const;
    double dNormal(double x) const;
};
#endif  // MATRIX_SURFACE_H_