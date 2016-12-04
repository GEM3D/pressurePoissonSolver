#include <vector>
using namespace std;
#include "Domain.h"

void solve(Domain &dom)
{
	vector<double> d_prime(dom.size());
	vector<double> c_prime(dom.size());

	// get some variables that are used for convenience
	const int    last_i = dom.size() - 1;
	const double h      = dom.spaceDelta();

	double left_boundary = 0.0;
	if (dom.hasLeftNbr()) {
		left_boundary = (dom.u.at(-1) + dom.u.at(0)) / 2;
	}
	double right_boundary = 0.0;
	if (dom.hasRightNbr()) {
		right_boundary = (dom.u.at(last_i) + dom.u.at(last_i + 1)) / 2;
	}

	// first value
	c_prime[0] = 1.0 / -3.0;
	d_prime[0] = (dom.uxx.at(0) * h * h - 2 * left_boundary) / (-3.0);
	// in between values
	for (int i = 1; i < last_i; i++) {
		c_prime[i] = 1.0 / (-2.0 - c_prime[i - 1]);
		d_prime[i] = (dom.uxx.at(i) * h * h - d_prime[i - 1]) / (-2.0 - c_prime[i - 1]);
	}
	// last value
	c_prime[last_i] = 1.0 / (-3.0 - c_prime[last_i - 1]);
	d_prime[last_i] = (dom.uxx.at(last_i) * h * h - 2 * right_boundary - d_prime[last_i - 1])
	                  / (-3.0 - c_prime[last_i - 1]);

	// back substitute
	dom.u.at(last_i) = d_prime[last_i];
	for (int i = last_i - 1; i >= 0; i--) {
		dom.u.at(i) = d_prime[i] - c_prime[i] * dom.u.at(i + 1);
	}
}
