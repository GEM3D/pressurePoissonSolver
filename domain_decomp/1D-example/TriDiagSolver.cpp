#include <iostream>
#include <valarray>
using namespace std;
#include "TriDiagSolver.h"

TriDiagSolver::TriDiagSolver(double left_boundary, double right_boundary)
{
	left_bnd  = left_boundary;
	right_bnd = right_boundary;
}

void TriDiagSolver::solve(Domain &dom)
{
	valarray<double> &u_xx = dom.u_xx;
	valarray<double> &u    = dom.u_curr;
	valarray<double>  d_prime(dom.size());
	valarray<double>  c_prime(dom.size());

	// get some variables that are used for convenience
	const int    last_i = dom.size() - 1;
	const double h      = dom.spaceDelta();

	// determine left boundary condition
	double left_boundary = left_bnd;
	if (dom.hasLeftNbr()) {
		left_boundary = dom.leftGamma();
	}

	// determine right boundary condition
	double right_boundary = right_bnd;
	if (dom.hasRightNbr()) {
		right_boundary = dom.rightGamma();
	}
	//cout << "On this domain using boundaries: " << left_boundary << ", " << right_boundary
	//     << "\n\n";

	// first value
	c_prime[0] = 1.0 / -3.0;
	d_prime[0] = (u_xx[0] * h * h - 2 * left_boundary) / -3.0;
	// in between values
	for (int i = 1; i < last_i; i++) {
		c_prime[i] = 1.0 / (-2.0 - c_prime[i - 1]);
		d_prime[i] = (u_xx[i] * h * h - d_prime[i - 1]) / (-2.0 - c_prime[i - 1]);
	}
	// last value
	c_prime[last_i] = 1.0 / (-3.0 - c_prime[last_i - 1]);
	d_prime[last_i] = (u_xx[last_i] * h * h - 2 * right_boundary - d_prime[last_i - 1])
	                  / (-3.0 - c_prime[last_i - 1]);

	// back substitute
	u[last_i] = d_prime[last_i];
	for (int i = last_i - 1; i >= 0; i--) {
		u[i] = d_prime[i] - c_prime[i] * u[i + 1];
	}
}
