#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace Eigen;
using namespace std;
#include "Domain.h"

typedef Matrix<Domain,-1,-1> DomainMatrix;
char *getCmdOption(char **begin, char **end, const std::string &option)
{
	char **itr = std::find(begin, end, option);
	if (itr != end && ++itr != end) {
		return *itr;
	}
	return 0;
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
	cout << "\t -s <file> \t save solution to file\n";
	cout << "\t -m \t print the matrix that was formed.\n";
	cout << "\t -h \t print this help message\n";
}

/**
 * @brief solve over all of the domains
 *
 * @param dmns the domains
 * @param gamma the interface values to use
 *
 * @return the difference between the interface values and gamma values.
 */
VectorXd solveWithInterfaceValues(DomainMatrix &dmns, VectorXd &gamma){
    //set the interface values on east and west
    int interface_size_ew = dmns(0,0).grid.rows();
	for (int j = 0; j < dmns.cols() - 1; j++) {
		for (int i = 0; i < dmns.rows(); i++) {
			int start_i = (j * dmns.rows() + i) * interface_size_ew;
			dmns(i, j).boundary_east     = gamma.block(start_i, 0, interface_size_ew, 1);
			dmns(i, j + 1).boundary_west = gamma.block(start_i, 0, interface_size_ew, 1);
		}
	}
    //set the interface values on north and south
    int interface_size_ns = dmns(0,0).grid.cols();
    int ns_start_i = dmns.rows()*(dmns.cols()-1)*interface_size_ew;
	for (int i = 0; i < dmns.rows()-1; i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			int start_i = (i * dmns.cols() + j) * interface_size_ns + ns_start_i;
			dmns(i, j).boundary_south = gamma.block(start_i, 0, interface_size_ns, 1).transpose();
			dmns(i + 1, j).boundary_north
			= gamma.block(start_i, 0, interface_size_ns, 1).transpose();
		}
	}

	//solve
	for (int i = 0; i < dmns.rows(); i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			dmns(i, j).solve();
		}
	}

    VectorXd diff(gamma.size());
	//get the difference on east-west interfaces
	for (int j = 0; j < dmns.cols() - 1; j++) {
		for (int i = 0; i < dmns.rows(); i++) {
			int start_i            = (j * dmns.rows() + i) * interface_size_ew;
			int left_dmns_last_col = dmns(i, j).grid.cols() - 1;
			diff.block(start_i, 0, interface_size_ew, 1)
			= dmns(i, j).u.col(left_dmns_last_col) + dmns(i, j + 1).u.col(0)
			  - 2 * gamma.block(start_i, 0, interface_size_ew, 1);
		}
	}
	//get the difference on north-south interfaces
	for (int i = 0; i < dmns.rows()-1; i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			int start_i           = (i * dmns.cols() + j) * interface_size_ns + ns_start_i;
			int top_dmns_last_row = dmns(i, j).grid.rows() - 1;
			diff.block(start_i, 0, interface_size_ns, 1)
			= dmns(i, j).u.row(top_dmns_last_row).transpose() + dmns(i+1, j).u.row(0).transpose()
			  - 2 * gamma.block(start_i, 0, interface_size_ns, 1);
		}
	}
	return diff;
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

// the functions that we are using
double ffun(double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); }
double gfun(double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); }

int main(int argc, char *argv[])
{
    // parse input
    string save_solution_file = "";
	bool   print_matrix       = false;
	if (argc < 2 || cmdOptionExists(argv, argv + argc, "-h")) {
		printHelp();
		return 1;
	}
	if (cmdOptionExists(argv, argv + argc, "-s")) {
		save_solution_file = getCmdOption(argv, argv + argc, "-s");
	}
	if (cmdOptionExists(argv, argv + argc, "-m")) {
		print_matrix = true;
	}

	int    m_x           = stoi(argv[1]);
	int    m_y           = m_x;
	int    num_domains_x = stoi(argv[2]);
	int    num_domains_y = num_domains_x;
	double x_start       = 0.0;
	double x_end         = 1.0;
	double h_x           = (x_end - x_start) / m_x;
	double h_y           = (x_end - x_start) / m_y;

	double error;

	// Populate the entire grid
	ArrayXf  x = ArrayXf::LinSpaced(m_x, x_start + h_x / 2.0, x_end - h_x / 2.0);
	ArrayXf  y = ArrayXf::LinSpaced(m_y, x_start + h_y / 2.0, x_end - h_y / 2.0);
	MatrixXd G(m_y, m_x);
	MatrixXd Exact(m_y, m_x);
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
		boundary_east(i) = gfun(x_start, y(i));
		boundary_west(i) = gfun(x_end, y(i));
	}

	RowVectorXd boundary_north(m_x);
	RowVectorXd boundary_south(m_x);
	for (int j = 0; j < m_x; j++) {
		boundary_north(j) = gfun(x(j), x_start);
		boundary_south(j) = gfun(x(j), x_end);
	}

	// Generate Domains
	DomainMatrix dmns(num_domains_y,num_domains_x);
	int            domain_width = m_x / num_domains_x;
	for (int i = 0; i < num_domains_y; i++) {
		for (int j = 0; j < num_domains_x; j++) {
			SparseMatrix<double> A        = generateMatrix(domain_width, domain_width, h_x, h_y);
			int                  start_i  = i * domain_width;
			int                  start_j  = j * domain_width;
			MatrixXd             sub_grid = G.block(start_i, start_j, domain_width, domain_width);
			dmns(i, j) = Domain(A, sub_grid, h_x, h_y);
		}
	}

	// set the outer boundary vectors on north and south
	for (int i = 0; i < num_domains_x; i++) {
		int start_j = i * m_x / num_domains_x;
		int length  = m_x / num_domains_x;
		dmns(0, i).boundary_north                 = boundary_north.block(0, start_j, 1, length);
		dmns(num_domains_y - 1, i).boundary_south = boundary_south.block(0, start_j, 1, length);
	}

	// set the outer boundary vectors on east and west
	for (int i = 0; i < num_domains_y; i++) {
		int start_i = i * m_y / num_domains_y;
		int length  = m_y / num_domains_y;
		dmns(i, 0).boundary_west                 = boundary_west.block(start_i,0, length, 1);
		dmns(i, num_domains_x - 1).boundary_east = boundary_east.block(start_i,0, length, 1);
	}

	// create gamma array
	VectorXd gamma            = VectorXd::Zero(m_y * (num_domains_x - 1)+m_x*(num_domains_y-1));
	double   condition_number = 0.0;
	if (num_domains_x > 1) {
		// get the b vector
		VectorXd b = solveWithInterfaceValues(dmns, gamma);

		// create a matrix
		MatrixXd A(gamma.size(), gamma.size());

		// get the columns of the matrix
		for (int i = 0; i < gamma.size(); i++) {
			gamma(i) = 1;
			A.col(i) = solveWithInterfaceValues(dmns, gamma) - b;
			gamma(i) = 0;
		}

		// solve for gamme values
		FullPivLU<MatrixXd> lu(A);
		VectorXd            tmp = lu.solve(b);
		gamma                   = -tmp;
		condition_number        = 1.0 / lu.rcond();

		if (print_matrix) {
			cout << "The matrix that was formed:\n";
			cout << A;
			cout << "\n\n";
		}
		if (save_solution_file != "") {
			// print out solution
			ofstream out_file(save_solution_file+".matrix");
			out_file.precision(20);
			out_file << scientific;
			out_file << A << "\n";
			out_file.close();
		}
	}

	// do one last solve
	solveWithInterfaceValues(dmns, gamma);

	// form the complete grid
	MatrixXd SOL(m_y, m_x);
	for (int i = 0; i < num_domains_y; i++) {
		for (int j = 0; j < num_domains_x; j++) {
			int start_i = i * domain_width;
			int start_j = j * domain_width;
			SOL.block(start_i, start_j, domain_width, domain_width) = dmns(i, j).u;
		}
	}

	// calculate error
	error = (SOL - Exact).norm() / Exact.norm();

	cout << scientific;
	cout.precision(13);
	cout << "Error: " << error << '\n';
	cout << defaultfloat;
	cout << "Condition Number: " << condition_number << '\n';

	if (save_solution_file != "") {
		// print out solution
		ofstream out_file(save_solution_file);
	    out_file.precision(20);
        out_file << scientific;
		out_file << SOL << "\n";
		out_file.close();
	}
}
