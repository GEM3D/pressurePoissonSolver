#include <cmath>
#include <iostream>
#include <limits>
#include <valarray>
#include <vector>
#define PI M_PI

using namespace std;
#include "Domain.h"
#include "TriDiagSolver.h"

// lapack function declaration
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);


double uxx_init(double x) { return -PI * PI * sin(PI * x); }
double exact_solution(double x) { return sin(PI * x); }
double error(vector<Domain> &dmns)
{
	double l2norm     = 0;
	double exact_norm = 0;
	for (Domain &d : dmns) {
		int    m       = d.size();
		double d_begin = d.domainBegin();
		double d_end   = d.domainEnd();
		for (int i = 0; i < m; i++) {
			double x     = d_begin + (i + 0.5) / m * (d_end - d_begin);
			double exact = exact_solution(x);
			double diff  = exact - d.u_curr[i];
			l2norm += diff * diff;
			exact_norm += exact * exact;
		}
	}
	return sqrt(l2norm) / sqrt(exact_norm);
}

/**
 * @brief solve over each of the domains
 *
 * @param tds the solver to use
 * @param dmns the domains to over
 * @param gammas the gamma values that are used
 *
 * @return
 */
valarray<double> solveOnAllDomains(TriDiagSolver &tds, vector<Domain> &dmns,
                                   valarray<double> &gammas)
{
	// solve over the domains
	for (Domain &d : dmns) {
		tds.solve(d);
	}

	// get the difference between the gamma value and computed solution at the interface
	valarray<double> z(gammas.size());
	for (size_t i = 0; i < gammas.size(); i++) {
		Domain &left_dmn  = dmns[i];
		Domain &right_dmn = dmns[i + 1];

		double left_val  = left_dmn.u_curr[left_dmn.u_curr.size() - 1];
		double right_val = right_dmn.u_curr[0];

		z[i] = left_val + right_val - 2 * gammas[i];
	}
	return z;
}

void printSolution(vector<Domain> &dmns)
{
	for (Domain &d : dmns) {
		for (double x : d.u_curr) {
			cout << x << "\t";
		}
	}
	cout << '\n';
}

int main(int argc, char *argv[])
{
	// set cout to print full precision
	// cout.precision(numeric_limits<double>::max_digits10);
	cout.precision(9);

	// create a solver with 0 for the boundary conditions
	TriDiagSolver  tds(0.0, 0.0);
	int            m           = stoi(argv[1]);
	int            num_domains = stoi(argv[2]);
	vector<Domain> dmns(num_domains);

	// create the domains
	for (int i = 0; i < num_domains; i++) {
		double x_start = (0.0 + i) / num_domains;
		double x_end   = (1.0 + i) / num_domains;
		dmns[i]        = Domain(x_start, x_end, m / num_domains, uxx_init);
	}

	// create an array to store the gamma values for each of the interfaces
	valarray<double> gammas(num_domains - 1);

	// set the gamma pointers
	if (num_domains > 1) {
		const int last_i        = num_domains - 1;
		dmns[0].right_gamma_ptr = &gammas[0];
		for (int i = 1; i < last_i; i++) {
			dmns[i].left_gamma_ptr  = &gammas[i - 1];
			dmns[i].right_gamma_ptr = &gammas[i];
		}
		dmns[last_i].left_gamma_ptr = &gammas[last_i - 1];
	}

	if (num_domains > 1) {
		// solve with gammas set to zero
		gammas             = 0;
		valarray<double> b = solveOnAllDomains(tds, dmns, gammas);

		cout << "b value(s):\n";
		for (double x : b) {
			cout << x << ' ';
		}
		cout << "\n\n";

		// build the A matrix
		int              n = gammas.size();
		valarray<double> A(n * n);
		for (int i = 0; i < n; i++) {
			gammas[i] = 1.0;
			A[slice(i * n, n, 1)] = solveOnAllDomains(tds, dmns, gammas) - b;
			gammas[i] = 0.0;
		}

		cout << "A matrix:\n";
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout << A[n * i + j] << '\t';
			}
			cout << '\n';
		}
		cout << "\n\n";

		// solve for the gamma values
		int           one = 1;
		valarray<int> ipiv(n);
		int           info;
		dgesv_(&n, &one, &A[0], &n, &ipiv[0], &b[0], &n, &info);

		// the solution gets stored in the b array
		gammas = -b;

		cout << "calculated gamma value(s):\n";
		for (double x : gammas) {
			cout << x << ' ';
		}
		cout << "\n\n";
	}

	/*
	 * get the final solution
	 */
	solveOnAllDomains(tds, dmns, gammas);

	// print out solution
	cout << "Final solution:\n";
	printSolution(dmns);
	cout << '\n';

	cerr << '\n';
	cerr << "error: " << scientific << error(dmns) << "\n";
}
