#include <vector>
using namespace std;
#include "TriDiagSolver.h"

TriDiagSolver::TriDiagSolver(double left_boundary, double right_boundary)
{
	left_bnd  = left_boundary;
	right_bnd = right_boundary;
}

void TriDiagSolver::solve(Domain &dom)
{
	vector<double> &uxx = dom.uxx();
	vector<double> &u   = dom.uCurr();
	vector<double>  d_prime(dom.size());
	vector<double>  c_prime(dom.size());

	// get some variables that are used for convenience
	const int    last_i = dom.size() - 1;
	const double h      = dom.spaceDelta();

	double left_boundary = left_bnd;
	if (dom.hasLeftNbr()) {
		left_boundary = (dom.leftNbr().uCurr().back() + u.at(0)) / 2;
	}
	double right_boundary = right_bnd;
	if (dom.hasRightNbr()) {
		right_boundary = (u.at(last_i) + dom.rightNbr().uPrev().front()) / 2;
	}

	// first value
	c_prime[0] = 1.0 / -3.0;
	d_prime[0] = (uxx.at(0) * h * h - 2 * left_boundary) / (-3.0);
	// in between values
	for (int i = 1; i < last_i; i++) {
		c_prime[i] = 1.0 / (-2.0 - c_prime[i - 1]);
		d_prime[i] = (uxx.at(i) * h * h - d_prime[i - 1]) / (-2.0 - c_prime[i - 1]);
	}
	// last value
	c_prime[last_i] = 1.0 / (-3.0 - c_prime[last_i - 1]);
	d_prime[last_i] = (uxx.at(last_i) * h * h - 2 * right_boundary - d_prime[last_i - 1])
	                  / (-3.0 - c_prime[last_i - 1]);

	// back substitute
	u.at(last_i) = d_prime[last_i];
	for (int i = last_i - 1; i >= 0; i--) {
		u.at(i) = d_prime[i] - c_prime[i] * u.at(i + 1);
	}
}
