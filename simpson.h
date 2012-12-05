#ifndef SIMPSON_H_
#define SIMPSON_H_

class Simpson{
  public:

  template <class F>
  static double Intergrate(double from, double to, int steps, const F& f){
    double  h = (to - from) / steps;
    double sum = f(from);

    for(int i = 1; i< steps; i++){
      if (i%2 == 0)
        sum += 2*f(from + h * i);
      else
        sum += 4*f(from + h * i);
    }

    sum += f(to);
    sum *= h/3;
    return sum;
  }
};

#endif  // SIMPSON_H_