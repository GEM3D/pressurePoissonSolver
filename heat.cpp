#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#define PI M_PI

using namespace std;
#include "Domain.h"
#include "TriDiagSolver.h"

double uxx_init(double x) { return -PI * PI * sin(PI * x); }
double exact_solution(double x) { return sin(PI * x); }
int main(int argc, char *argv[])
{
	// set cout to print full precision
	cout.precision(numeric_limits<double>::max_digits10);

	// create a solver with 0 for the boundary conditions
	TriDiagSolver    tds(0.0, 0.0);
	int              m           = stoi(argv[1]);
	int              num_domains = stoi(argv[2]);
	vector<Domain *> dmns(num_domains);

	// create the domains
	for (int i = 0; i < num_domains; i++) {
		double x_start = (0.0 + i) / num_domains;
		double x_end   = (1.0 + i) / num_domains;
		dmns[i] = new Domain(x_start, x_end, m / num_domains, uxx_init);
	}

	// set neighbors
	if (num_domains > 1) {
		const int last_i = num_domains - 1;
		dmns[0]->setRightNbr(*dmns[1]);
		for (int i = 1; i < last_i; i++) {
			dmns[i]->setLeftNbr(*dmns[i - 1]);
			dmns[i]->setRightNbr(*dmns[i + 1]);
		}
		dmns[last_i]->setLeftNbr(*dmns[last_i - 1]);
	}

	// get l2norm of uxx
	double uxx_l2norm = 0;
	for (Domain *d_ptr : dmns) {
		for (double x : d_ptr->uxx()) {
			uxx_l2norm += x * x;
		}
	}
	uxx_l2norm = sqrt(uxx_l2norm);

	// start solving
	double l2norm   = 0;
	int    num_iter = 0;
	do {
		num_iter++;

		// solve over each domain
		for (Domain *d_ptr : dmns) {
			tds.solve(*d_ptr);
		}

		// print out solution
		for (Domain *d_ptr : dmns) {
			for (double x : d_ptr->uCurr()) {
				cout << x << "\t";
			}
		}
		cout << '\n';

		// calculate l2norm
		if (num_domains > 1) {
			l2norm = 0;
			for (Domain *d_ptr : dmns) {
				double tmp = d_ptr->l2norm();
				l2norm += tmp * tmp;
				// go ahead and swap the vectors
				d_ptr->swapCurrPrev();
			}
			l2norm = sqrt(l2norm);
		}
	} while (l2norm > uxx_l2norm * 10e-10);

	// delete the domains
	for (Domain *d_ptr : dmns) {
		delete d_ptr;
	}
	cerr << '\n';
	cerr << "number of iterations: " << num_iter << "\n";
}
