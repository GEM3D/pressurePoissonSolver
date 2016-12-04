#include <cmath>
#include <iostream>
#include <vector>
#define PI M_PI

using namespace std;
#include "Grid.h"
#include "TriDiagSolver.h"

double f(double x) { return -PI * PI * sin(PI * x); }
int             main()
{
	int  m = 10;
	Grid left_grid(0.0, 0.5, m / 2, f);
	Grid right_grid(0.5, 1.0, m / 2, f);
	// set neighbors
	left_grid.setRightNbr(right_grid);
	right_grid.setLeftNbr(left_grid);
	for (int i = 0; i < left_grid.size(); i++) {
		cout << left_grid.uxx.at(i) << "  ";
	}
	for (int i = 0; i < right_grid.size(); i++) {
		cout << right_grid.uxx.at(i) << "  ";
	}
	cout << '\n';
	for (int i = 0; i < 100; i++) {
		solve(left_grid);
		solve(right_grid);
		for (int i = 0; i < left_grid.size(); i++) {
			cout << left_grid.u.at(i) << "  ";
		}
		for (int i = 0; i < right_grid.size(); i++) {
			cout << right_grid.u.at(i) << "  ";
		}
		cout << '\n';
	}
	cout << '\n';
	cout << "The soluiton on a single grid:\n";
	Grid single_grid(0.0, 1.0, m, f);
	solve(single_grid);
	for (int i = 0; i < single_grid.size(); i++) {
		cout << single_grid.u.at(i) << "  ";
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
