#ifndef SIMPSON_H_
#define SIMPSON_H_

#include "utils.h"
namespace{

template<typename F>
class SubSum {
public:
  SubSum(const F& f, double offset, double scale)
    : f_(f), offset_(offset), scale_(scale) {
  }

  double operator () (int i) const {
    return f_(offset_ + scale_ * i);
  }

private:
  const F& f_;
  double offset_;
  double scale_;
};
}

class Simpson{
  public:

  template <class F>
  static double Integrate(double from, double to, int steps, const F& f){
    double  h = (to - from) / steps;

    SubSum<F> even_sum(f, from + 2 * h, 2 * h);
    SubSum<F> odd_sum(f, from + h, 2 * h);

    double sum1 = utils::KahanSum((steps - 1) / 2, even_sum);
    double sum2 = utils::KahanSum(steps / 2, odd_sum);

    double sum = f(from) + 2*sum1 + 4*sum2 + f(to);
    sum *= h/3;
    return sum;
  }
};

#endif  // SIMPSON_H_