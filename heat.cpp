#include <cmath>
#include <iostream>
#include <vector>
#define PI M_PI

using namespace std;
#include "Domain.h"
#include "TriDiagSolver.h"

double uxx_init(double x) { return -PI * PI * sin(PI * x); }
int             main()
{
	TriDiagSolver tds(0.0, 0.0);
	int           m           = 10;
	int           num_domains = 2;
	Domain        left_dom(0.0, 0.5, m / 2, uxx_init);
	Domain        right_dom(0.5, 1.0, m / 2, uxx_init);
	// set neighbors
	left_dom.setRightNbr(right_dom);
	right_dom.setLeftNbr(left_dom);
	cout << '\n';
	double l2norm = 0;
	do {
		tds.solve(left_dom);
		tds.solve(right_dom);
		for (double x : left_dom.uCurr()) {
			cout << x << "  ";
		}
		for (double x : right_dom.uCurr()) {
			cout << x << "  ";
		}
		cout << '\n';
		left_dom.swapCurrPrev();
		right_dom.swapCurrPrev();
		double tmp = left_dom.l2norm();
		l2norm     = tmp * tmp;
		tmp        = right_dom.l2norm();
		l2norm += tmp * tmp;
		l2norm = sqrt(l2norm);
	} while (l2norm > 10e-10);
	cout << '\n';
}
