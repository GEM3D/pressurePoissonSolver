#include <cmath>
#include <iostream>
#include <vector>
#define PI M_PI

using namespace std;
#include "Domain.h"
#include "TriDiagSolver.h"

double f(double x) { return -PI * PI * sin(PI * x); }
int             main()
{
	int    m = 10;
	Domain left_dom(0.0, 0.5, m / 2, f);
	Domain right_dom(0.5, 1.0, m / 2, f);
	// set neighbors
	left_dom.setRightNbr(right_dom);
	right_dom.setLeftNbr(left_dom);
	for (double x : left_dom.getGrid("u_xx")) {
		cout << x << "  ";
	}
	for (double x : right_dom.getGrid("u_xx")) {
		cout << x << "  ";
	}
	cout << '\n';
	for (int i = 0; i < 100; i++) {
		solve(left_dom);
		solve(right_dom);
		for (double x : left_dom.getGrid("u")) {
			cout << x << "  ";
		}
		for (double x : right_dom.getGrid("u")) {
			cout << x << "  ";
		}
		cout << '\n';
	}
	cout << '\n';
	cout << "The soluiton on a single dom:\n";
	Domain single_dom(0.0, 1.0, m, f);
	solve(single_dom);
	for (double x : single_dom.getGrid("u")) {
		cout << x << "  ";
	}
	cout << '\n' << '\n';

	// what the exct solution should be
	cout << "The exact solution:\n";
	for (int i = 0; i < m; i++) {
		double x      = (0.5 + i) / m;
		double result = sin(PI * x);
		cout << result << "  ";
	}
	cout << '\n';
}
