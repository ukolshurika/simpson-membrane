#ifndef MATRIX_SURFACE_H_
#define MATRIX_SURFACE_H_

class MatrixSurface{
  public:
    double operator () (double x) const;
    double derivative (double x) const;
    double alphaConstrained () const;
    double right_zero () const;
};
#endif  // MATRIX_SURFACE_H_