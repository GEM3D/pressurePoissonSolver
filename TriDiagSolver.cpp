#include <vector>
using namespace std;
#include "Grid.h"

void solve(Grid &grid)
{
	vector<double> d_prime(grid.size());
	vector<double> c_prime(grid.size());

	// get some variables that are used for convenience
	const int    last_i = grid.size() - 1;
	const double h      = grid.spaceDelta();

	double left_boundary = 0.0;
	if (grid.hasLeftNbr()) {
		left_boundary = (grid.u.at(-1) + grid.u.at(0)) / 2;
	}
	double right_boundary = 0.0;
	if (grid.hasRightNbr()) {
		right_boundary = (grid.u.at(last_i) + grid.u.at(last_i + 1)) / 2;
	}

	// first value
	c_prime[0] = 1.0 / -3.0;
	d_prime[0] = (grid.uxx.at(0) * h * h - 2 * left_boundary) / (-3.0);
	// in between values
	for (int i = 1; i < last_i; i++) {
		c_prime[i] = 1.0 / (-2.0 - c_prime[i - 1]);
		d_prime[i] = (grid.uxx.at(i) * h * h - d_prime[i - 1]) / (-2.0 - c_prime[i - 1]);
	}
	// last value
	c_prime[last_i] = 1.0 / (-3.0 - c_prime[last_i - 1]);
	d_prime[last_i] = (grid.uxx.at(last_i) * h * h - 2 * right_boundary - d_prime[last_i - 1])
	                  / (-3.0 - c_prime[last_i - 1]);

	// back substitute
	grid.u.at(last_i) = d_prime[last_i];
	for (int i = last_i - 1; i >= 0; i--) {
		grid.u.at(i) = d_prime[i] - c_prime[i] * grid.u.at(i + 1);
	}
}
