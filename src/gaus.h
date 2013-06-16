#ifndef SIMPSON_H_
#define SIMPSON_H_

#include <cmath>

class Gaus{
  public:

  template <class F>
  static double Integrate(double from, double to, const F& f){
    double summand = (from+to)/2.0, sum, multiplier = (to-from)/2.0;
    sum =  128/255.0*f(summand)+
           (322+13*sqrt(70))/900.0*f(summand+1/3.0*sqrt(5-2*sqrt(10/7.0)))*multiplier +
           (322+13*sqrt(70))/900.0*f(summand-1/3.0*sqrt(5-2*sqrt(10/7.0)))*multiplier +
           (322-13*sqrt(70))/900.0*f(summand+1/3.0*sqrt(5+2*sqrt(10/7.0)))*multiplier +
           (322-13*sqrt(70))/900.0*f(summand-1/3.0*sqrt(5+2*sqrt(10/7.0)))*multiplier;

    return sum*(to-from)/2.0;
  }
};

#endif  // SIMPSON_H_