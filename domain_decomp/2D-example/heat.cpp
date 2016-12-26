#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;
#include "Domain.h"
double ffun(double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); }
double gfun(double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); }
void solveWithInterfaceValues(vector<Domain> &dmns){
    return;
}
SparseMatrix<double> generateMatrix(int m_x, int m_y, double h_x, double h_y){
	// crate an identity matricies
	SparseMatrix<double> eye_x(m_x, m_x);
	eye_x.setIdentity();
	SparseMatrix<double> eye_y(m_y, m_y);
	eye_y.setIdentity();

	// create diagonal for x values
	SparseMatrix<double> Diag_x(m_x, m_x);
	Diag_x.insert(0, 0) = -3 / (h_x * h_x);
	for (int i = 1; i < m_x - 1; i++) {
		Diag_x.insert(i, i)     = -2 / (h_x * h_x);
		Diag_x.insert(i - 1, i) = 1 / (h_x * h_x);
		Diag_x.insert(i, i - 1) = 1 / (h_x * h_x);
	}
	if (m_x > 1) {
		int i = m_x - 1;
		Diag_x.insert(i, i)     = -3 / (h_x * h_x);
		Diag_x.insert(i - 1, i) = 1 / (h_x * h_x);
		Diag_x.insert(i, i - 1) = 1 / (h_x * h_x);
	}

	// create diagonal for y values
	SparseMatrix<double> Diag_y(m_y, m_y);
	Diag_y.insert(0, 0) = -3 / (h_y * h_y);
	for (int i = 1; i < m_y - 1; i++) {
		Diag_y.insert(i, i)     = -2 / (h_y * h_y);
		Diag_y.insert(i - 1, i) = 1 / (h_y * h_y);
		Diag_y.insert(i, i - 1) = 1 / (h_y * h_y);
	}
	if (m_y > 1) {
		int i = m_y - 1;
		Diag_y.insert(i, i)     = -3 / (h_y * h_y);
		Diag_y.insert(i - 1, i) = 1 / (h_y * h_y);
		Diag_y.insert(i, i - 1) = 1 / (h_y * h_y);
	}
	// form the matrix
	SparseMatrix<double> A = kroneckerProduct(Diag_x, eye_y) + kroneckerProduct(eye_x, Diag_y);
    return A;
}
int main(int argc, char *argv[])
{
	int    m_x         = 512;
	int    m_y         = 512;
	int    num_domains = 1;
	double x_start = 0.0;
	double x_end   = 1.0;
	double h_x       = (x_end - x_start) / m_x;
	double h_y       = (x_end - x_start) / m_y;

	double error;
	ArrayXf  x = ArrayXf::LinSpaced(m_x, x_start + h_x / 2.0, x_end - h_x / 2.0);
	ArrayXf  y = ArrayXf::LinSpaced(m_y, x_start + h_y / 2.0, x_end - h_y / 2.0);
	MatrixXd G(m_y, m_x);
	MatrixXd Exact(m_y, m_x);
	// Generate the entire grid
	for (int i = 0; i < m_y; i++) {
		for (int j = 0; j < m_x; j++) {
			G(i, j)     = ffun(x(j), y(i));
			Exact(i, j) = gfun(x(j), y(i));
		}
	}

	// generate boundary vectors
    VectorXd boundary_east(m_y);
    VectorXd boundary_west(m_y);
	for (int i = 0; i < m_y; i++) {
		boundary_east(i)       = gfun(x_start, y(i));
		boundary_west(i)       = gfun(x_end, y(i));
	}
    RowVectorXd boundary_north(m_x);
    RowVectorXd boundary_south(m_x);
	for (int j = 0; j < m_x; j++) {
		boundary_north(j) = gfun(x(j), x_start);
		boundary_south(j) = gfun(x(j), x_end);
	}

    SparseMatrix<double> A = generateMatrix(m_x,m_y,h_x,h_y);
	Domain d = Domain(A, G, h_x, h_y);
    // set the boundary vectors
    d.boundary_north = boundary_north;
    d.boundary_south = boundary_south;
    d.boundary_east = boundary_east;
    d.boundary_west = boundary_west;
	d.solve();

	// cout << A << "\n\n";
	// cout << G << "\n\n";
	// cout << Exact << "\n\n";

	// SimplicialLDLT<SparseMatrix<double>> solver;
	// solver.compute(A);
	// VectorXd      xSol = solver.solve(f);
	MatrixXd SOL = d.u;
	error        = (SOL - Exact).norm() / Exact.norm();

	// cout << SOL << "\n\n";
	cout << "Error: " << error << '\n';
}
