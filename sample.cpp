#include <cmath>
#include <iostream>
using namespace std;

struct Square {
	double operator () (double x) const {
		return x * x;
	}
};

struct SquareRoot {
	double operator () (double x) const {
		return sqrt(x);
	}
};

template <class F>
void Calculate(const F& f) {
	std::cout << f(10.0) << std::endl;
}


int main() {
	Calculate(Square());
	Calculate(SquareRoot());
	return 0;
}
