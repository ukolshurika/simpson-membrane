#ifndef FREE_DEFORMATION_H_
#define FREE_DEFORMATION_H_


class FreeDeformation{
  public:
  FreeDeformation(double h0, double q, double n);
  double operator () (double alpha) const;
  double h(double alpha);

  double h0_;
  double q_;
  double n_;
};

#endif  // FREE_DEFORMATION_H