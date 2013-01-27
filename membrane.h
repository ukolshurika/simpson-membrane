#ifndef MEMBRANE_H_
#define MEMBRANE_H_

#include <map>
#include <vector>
#include <utility>
#include <cmath>

#include "matrix_surface.h"


class Membrane{
  public:
    Membrane(double h0, double q, double n, double sigma_b, double a,  double epsilon, int simpsonStep, int steps);
    void EasyIntegrate(int tid, std::vector<std::pair<double, double>>* v);
    void IntegrateConstrained(std::vector<std::pair<double, double>>* v);
    double operator () (double alpha) const;
    void IntegrateForAnimation(int steps);
    void OutputResult();
    void FreeStep(int steps);
    void ConstrainedStep(int steps);
  // private:
    std::map<double, double> times_free_;
    std::map<double, double> times_constrained_;

    double h0_;
    double q_;
    double n_;
    double sigma_b_;
    double a_;
    double epsilon_;
    int dstep_;
    double da_;
    double dx_;
    int simpsonStep_;
    int steps_;
    int num_threads_;
    MatrixSurface m_surface_;

    double ValueAsLine(double time, std::map<double, double>::iterator point1, std::map<double, double>::iterator point2);
    void CorrectTimes();
    void AverageDt(double dt_average);
    double MeanValueDt();
    double MeanValueDa();
};
#endif  // MEMBRANE_H_