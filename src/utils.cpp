#include "utils.h"

bool utils::eql(double a, double b){
	return abs(a-b) < kEpsilon;
};

bool utils::IsNan(double a){
	return a != a
};
