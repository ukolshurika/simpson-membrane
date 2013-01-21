#ifndef MATRIX_SURFACE_H_
#define MATRIX_SURFACE_H_

#define DELTA 0.001
class MatrixSurface{
  public:
    double operator () (double x) const;
    double derivative (double x) const;
    double alphaConstrained () const;
    double right_zero () const;
    double 
};
#endif  // MATRIX_SURFACE_H_