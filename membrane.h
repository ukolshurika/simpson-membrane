#ifndef MEMBRANE_H_
#define MEMBRANE_H_

#include <map>

class Membrane{
  public:
    Membrane(double h0, double q, double n, double epsilon, int simpsonStep, int steps);
    void EasyIntegrate(int steps);
    double operator () (double alpha) const;
    void IntegrateForAnimation(int steps);
    void OutputResult();
  private:
    std::map<double, double> times_;
    std::map<double, double> da_;

    double h0_;
    double q_;
    double n_;
    double epsilon_;
    double dstep_
    double da_;
    int steps_;
    int simpsonStep_;
    int num_theads;
    double ValueAsLine(double time, std::map<double, double>::iterator point1, std::map<double, double>::iterator point2);
    void CorrectTimes();
    void AverageDt(double dt_average);
    double MeanValueDt();
    double MeanValueDa();
};
#endif  // MEMBRANE_H_