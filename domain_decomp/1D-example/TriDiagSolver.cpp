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

	// determine left boundary conditions
	double left_boundary  = 2 * left_bnd;
	double top_left_coeff = -3.0;
	if (dom.hasLeftNbr()) {
		top_left_coeff             = -2.0;
		valarray<double> &left_sol = dom.leftNbr().u_curr;
		left_boundary              = left_sol[left_sol.size() - 1];
	}

	// determine right boundary conditions
	double right_boundary     = 2 * right_bnd;
	double bottom_right_coeff = -3.0;
	if (dom.hasRightNbr()) {
		bottom_right_coeff          = -2.0;
		valarray<double> &right_sol = dom.rightNbr().u_prev;
		right_boundary              = right_sol[0];
	}

	// first value
	c_prime[0] = 1.0 / top_left_coeff;
	d_prime[0] = (u_xx[0] * h * h - left_boundary) / top_left_coeff;
	// in between values
	for (int i = 1; i < last_i; i++) {
		c_prime[i] = 1.0 / (-2.0 - c_prime[i - 1]);
		d_prime[i] = (u_xx[i] * h * h - d_prime[i - 1]) / (-2.0 - c_prime[i - 1]);
	}
	// last value
	c_prime[last_i] = 1.0 / (bottom_right_coeff - c_prime[last_i - 1]);
	d_prime[last_i] = (u_xx[last_i] * h * h - right_boundary - d_prime[last_i - 1])
	                  / (bottom_right_coeff - c_prime[last_i - 1]);

	// back substitute
	u[last_i] = d_prime[last_i];
	for (int i = last_i - 1; i >= 0; i--) {
		u[i] = d_prime[i] - c_prime[i] * u[i + 1];
	}
}
