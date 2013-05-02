#ifndef MEMBRANE_H_
#define MEMBRANE_H_

#include <vector>
#include <utility>
#include <string>

#include "matrix.h"

class Membrane{
  public:
  Membrane(double q, double h0, double n);

  void free(int steps);
  void constrained(int steps);
  void averageAllDt();
  
  std::string free_data_to_draw(double alpha);
  std::string constrained_x_data_to_draw(double x);
  std::string constrained_y_data_to_draw(double x);


  std::vector<std::pair<double, double>> t_free_;
  std::vector<std::pair<double, double>> t_constrained_;
  std::vector<std::pair<double, double>> t_constrained_y_;
  double q_;
  double n_;


  /*величины толщины после первого и второго этапа*/
  double h1_;
  double h2_;

  // private:
  double h0_;


  /*величина угла, после упругого деформирования*/
  double alpha1_;
  /*величина угла, при котором происходит касание*/
  double alpha2_;
};



#endif