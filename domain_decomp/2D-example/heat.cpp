#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace Eigen;
using namespace std;
#include "Domain.h"

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
VectorXd solveWithInterfaceValues(vector<Domain> &dmns, VectorXd &gamma){
    //set the interface values
    int interface_size = dmns[0].grid.rows();
	for (size_t i = 0; i < dmns.size()-1; i++) {
        int start_i = i*interface_size;
        dmns[i].boundary_east = gamma.block(start_i,0,interface_size,1);
        dmns[i+1].boundary_west = gamma.block(start_i,0,interface_size,1);
    }
    
    //solve
    for(Domain& d: dmns){
        d.solve();
    }

    //get the difference
    VectorXd diff(gamma.size());
	for (size_t i = 0; i < dmns.size()-1; i++) {
        int start_i = i*interface_size;
        int left_dmns_last_col = dmns[i].grid.cols()-1;
		diff.block(start_i, 0, interface_size, 1)
		= dmns[i].u.col(left_dmns_last_col) + dmns[i + 1].u.col(0)
		  - 2 * gamma.block(start_i, 0, interface_size, 1);
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

	int    m_x         = stoi(argv[1]);
	int    m_y         = m_x;
	int    num_domains = stoi(argv[2]);
	double x_start     = 0.0;
	double x_end       = 1.0;
	double h_x         = (x_end - x_start) / m_x;
	double h_y         = (x_end - x_start) / m_y;

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
	vector<Domain> dmns(num_domains);
	int            domain_width = m_x / num_domains;
	for (int i = 0; i < num_domains; i++) {
		SparseMatrix<double> A        = generateMatrix(domain_width, m_y, h_x, h_y);
		int                  start_j  = i * m_x / num_domains;
		MatrixXd             sub_grid = G.block(0, start_j, m_y, domain_width);
		dmns[i]                       = Domain(A, sub_grid, h_x, h_y);
	}

	// set the outer boundary vectors
	for (int i = 0; i < num_domains; i++) {
		if (i == 0) {
			dmns[0].boundary_west = boundary_west;
		}
		if (i == num_domains - 1) {
			dmns[num_domains - 1].boundary_east = boundary_east;
		}
		int start_j            = i * m_x / num_domains;
		int length             = m_x / num_domains;
		dmns[i].boundary_north = boundary_north.block(0, start_j, 1, length);
		dmns[i].boundary_south = boundary_south.block(0, start_j, 1, length);
	}

	// create gamma array
	VectorXd gamma            = VectorXd::Zero(m_y * (num_domains - 1));
	double   condition_number = 0.0;
	if (num_domains > 1) {
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
	}

	// do one last solve
	solveWithInterfaceValues(dmns, gamma);

	// form the complete grid
	MatrixXd SOL(m_y, m_x);
	for (int i = 0; i < num_domains; i++) {
		int start_j = i * m_x / num_domains;
		SOL.block(0, start_j, m_y, domain_width) = dmns[i].u;
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
		out_file << SOL << "\n";
		out_file.close();
	}
}
