#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;
double ffun(double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); }
double gfun(double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); }
int main(int argc, char *argv[])
{
	int    m       = stoi(argv[1]);
	double x_start = 0.0;
	double x_end   = 1.0;
	double h       = (x_end - x_start) / (m + 1);

	ArrayXf  x = ArrayXf::LinSpaced(m + 2, x_start, x_end);
	ArrayXf  y = ArrayXf::LinSpaced(m + 2, x_start, x_end);
	MatrixXd G(m, m);
	MatrixXd Exact(m, m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			G(i, j)     = ffun(x(i + 1), y(j + 1));
			Exact(i, j) = gfun(x(i + 1), y(j + 1));
		}
	}
	for (int i = 0; i < m; i++) {
		G(i, 0)     = G(i, 0) - gfun(x(i + 1), y(0)) / (h * h);
		G(i, m - 1) = G(i, m - 1) - gfun(x(i + 1), y(m + 1)) / (h * h);
	}
	for (int j = 0; j < m; j++) {
		G(0, j)     = G(0, j) - gfun(x(0), y(j + 1)) / (h * h);
		G(m - 1, j) = G(m - 1, j) - gfun(x(m + 1), y(j + 1)) / (h * h);
	}

	Map<VectorXd> f(G.data(), G.size());

	// crate an identity matrix
	SparseMatrix<double> eye(m, m);
	eye.setIdentity();

	SparseMatrix<double> D2(m, m);
	D2.insert(0, 0) = -2 / (h * h);
	for (int i = 1; i < m; i++) {
		D2.insert(i, i)     = -2 / (h * h);
		D2.insert(i - 1, i) = 1 / (h * h);
		D2.insert(i, i - 1) = 1 / (h * h);
	}
	SparseMatrix<double> A = kroneckerProduct(D2, eye) + kroneckerProduct(eye, D2);

	// cout << A << "\n\n";
	// cout << G << "\n\n";
	// cout << Exact << "\n\n";

	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(A);
	VectorXd      xSol = solver.solve(f);
	Map<MatrixXd> SOL(xSol.data(), m, m);

	// cout << SOL << "\n\n";
	cout << "Error: " << (SOL - Exact).norm() / Exact.norm() << '\n';
}
