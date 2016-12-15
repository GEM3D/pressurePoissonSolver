#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <valarray>
#include <vector>
#define PI M_PI

using namespace std;
#include "Domain.h"
#include "TriDiagSolver.h"

// lapack function signatures
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb,
                       int *info);
extern "C" void dgecon_(char *norm, int *n, double *A, int *lda, double *Anorm, double *rcond,
                        double *work, int *iwork, int *info);

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
 *
 * @return the difference between the gamma values and the computed value at the domain
 */
valarray<double> solveOnAllDomains(TriDiagSolver &tds, vector<Domain> &dmns)
{
	// solve over the domains
	for (Domain &d : dmns) {
		tds.solve(d);
	}

	// get the difference between the gamma value and computed solution at the interface
	valarray<double> z(dmns.size() - 1);
	for (size_t i = 0; i < z.size(); i++) {
		Domain &left_dmn  = dmns[i];
		Domain &right_dmn = dmns[i + 1];
		double  gamma     = left_dmn.rightGamma();

		double left_val  = left_dmn.u_curr[left_dmn.u_curr.size() - 1];
		double right_val = right_dmn.u_curr[0];

		z[i] = left_val + right_val - 2 * gamma;
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

bool cmdOptionExists(char **begin, char **end, const string &option)
{
	return find(begin, end, option) != end;
}

void printHelp()
{
	cout << "Usage:\n";
	cout << "./heat <num_cells> <numdomains> [options]\n\n";
	cout << "Options can be:\n";
	cout << "\t -s \t print solution\n";
	cout << "\t -m \t print the matrix that was formed.\n";
	cout << "\t -h \t print this help message\n";
}

int main(int argc, char *argv[])
{
	bool print_solution = false;
	bool print_matrix   = false;
	if (argc < 2 || cmdOptionExists(argv, argv + argc, "-h")) {
		printHelp();
		return 1;
	}
	if (cmdOptionExists(argv, argv + argc, "-s")) {
		print_solution = true;
	}
	if (cmdOptionExists(argv, argv + argc, "-m")) {
		print_matrix = true;
	}

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

	double condition_number;
	if (num_domains > 1) {
		// solve with gammas set to zero
		gammas             = 0;
		valarray<double> b = solveOnAllDomains(tds, dmns);

		if (print_matrix) {
			cout << "b value(s):\n";
			for (double x : b) {
				cout << x << ' ';
			}
			cout << "\n\n";
		}

		// build the A matrix
		int              n = gammas.size();
		valarray<double> A(n * n);
		for (int i = 0; i < n; i++) {
			gammas[i] = 1.0;
			A[slice(i * n, n, 1)] = solveOnAllDomains(tds, dmns) - b;
			gammas[i] = 0.0;
		}

		if (print_matrix) {
			cout << "A matrix:\n";
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					cout << A[n * i + j] << '\t';
				}
				cout << '\n';
			}
			cout << "\n\n";
		}

		// get the condition number of A
		double rcond;

		{
			// I'm putting these in their own scope so I don't clutter up the namespace with these
			// variables
			char             norm = '1';
			valarray<double> col_sum(n);
			for (int i = 0; i < n; i++) {
				valarray<double> col = A[slice(i * n, n, 1)];
				col_sum[i]           = abs(col).sum();
			}
			double           A_norm = col_sum.max();
			valarray<double> work(4 * n);
			valarray<int>    iwork(n);
			int              info;
			dgecon_(&norm, &n, &A[0], &n, &A_norm, &rcond, &work[0], &iwork[0], &info);
			assert(info == 0);
		}
		condition_number = 1.0 / rcond;

		// solve for the gamma values
		{
			int           one = 1;
			valarray<int> ipiv(n);
			int           info;
			dgesv_(&n, &one, &A[0], &n, &ipiv[0], &b[0], &n, &info);
			assert(info == 0);
		}

		// the solution gets stored in the b array
		for (int i = 0; i < n; i++) {
			gammas[i] = -b[i];
		}

		if (print_matrix) {
			cout << "calculated gamma value(s):\n";
			for (double x : gammas) {
				cout << x << ' ';
			}
			cout << "\n\n";
		}
	}

	/*
	 * get the final solution
	 */
	solveOnAllDomains(tds, dmns);

	if (print_solution) {
		// print out solution
		cout << "Final solution:\n";
		printSolution(dmns);
		cout << '\n';
	}

	cout << "error: " << scientific << error(dmns) << "\n" << defaultfloat;
	if (num_domains > 1) {
		cout << "condition number of A: " << condition_number << "\n";
	}
}
