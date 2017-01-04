#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/KroneckerProduct>
using namespace Eigen;
using namespace std;
#include "Domain.h"
#include "args.h"

typedef Matrix<Domain, -1, -1> DomainMatrix;

class FunctionWrapper;
using Eigen::SparseMatrix;
namespace Eigen
{
namespace internal
{
// FunctionWrapper looks-like a SparseMatrix, so let's inherits its traits:
template <>
struct traits<FunctionWrapper> : public Eigen::internal::traits<Eigen::SparseMatrix<double> > {
};
}
}
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class FunctionWrapper : public Eigen::EigenBase<FunctionWrapper>
{
	private:
	const SparseMatrix<double> *mp_mat;
	Index                       my_rows;

	public:
	DomainMatrix *dmns;
	VectorXd      b;
	// Required typedefs, constants, and method:
	typedef double Scalar;
	typedef double RealScalar;
	typedef int    StorageIndex;
	enum {
		ColsAtCompileTime    = Eigen::Dynamic,
		MaxColsAtCompileTime = Eigen::Dynamic,
		IsRowMajor           = false
	};
	Index rows() const { return my_rows; }
	Index cols() const { return my_rows; }
	template <typename Rhs>
	Eigen::Product<FunctionWrapper, Rhs, Eigen::AliasFreeProduct>
	operator*(const Eigen::MatrixBase<Rhs> &x) const
	{
		return Eigen::Product<FunctionWrapper, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}
	// Custom API:
	FunctionWrapper(DomainMatrix *dmns, int size, VectorXd b)
	{
		this->dmns = dmns;
		this->b    = b;
		my_rows    = size;
	}

	void attachMyMatrix(const SparseMatrix<double> &mat) { mp_mat = &mat; }
	const SparseMatrix<double>                      my_matrix() const { return *mp_mat; }
};
VectorXd solveWithInterfaceValues(DomainMatrix &dmns, VectorXd &gamma);
// Implementation of FunctionWrapper * Eigen::DenseVector though a specialization of
// internal::generic_product_impl:
namespace Eigen
{
namespace internal
{
template <typename Rhs>
struct generic_product_impl<FunctionWrapper, Rhs, SparseShape, DenseShape,
                            GemvProduct> // GEMV stands for matrix-vector
: generic_product_impl_base<FunctionWrapper, Rhs, generic_product_impl<FunctionWrapper, Rhs> > {
	typedef typename Product<FunctionWrapper, Rhs>::Scalar Scalar;
	template <typename Dest>
	static void scaleAndAddTo(Dest &dst, const FunctionWrapper &lhs, const Rhs &rhs,
	                          const Scalar &alpha)
	{
		// This method should implement "dst += alpha * lhs * rhs" inplace,
		// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
		assert(alpha == Scalar(1) && "scaling is not implemented");
		// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
		// but let's do something fancier (and less efficient):
		VectorXd rhscopy = rhs;
		VectorXd result  = solveWithInterfaceValues(*lhs.dmns, rhscopy) - lhs.b;
		dst += alpha * result;
	}
};
}
}

/**
 * @brief solve over all of the domains
 *
 * @param dmns the domains
 * @param gamma the interface values to use
 *
 * @return the difference between the interface values and gamma values.
 */
VectorXd solveWithInterfaceValues(DomainMatrix &dmns, VectorXd &gamma)
{
	// set the interface values on east and west
	int interface_size_ew = dmns(0, 0).grid.rows();
	for (int j = 0; j < dmns.cols() - 1; j++) {
		for (int i = 0; i < dmns.rows(); i++) {
			int start_i = (j * dmns.rows() + i) * interface_size_ew;
			dmns(i, j).boundary_east     = gamma.block(start_i, 0, interface_size_ew, 1);
			dmns(i, j + 1).boundary_west = gamma.block(start_i, 0, interface_size_ew, 1);
		}
	}
	// set the interface values on north and south
	int interface_size_ns = dmns(0, 0).grid.cols();
	int ns_start_i        = dmns.rows() * (dmns.cols() - 1) * interface_size_ew;
	for (int i = 0; i < dmns.rows() - 1; i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			int start_i = (i * dmns.cols() + j) * interface_size_ns + ns_start_i;
			dmns(i, j).boundary_south = gamma.block(start_i, 0, interface_size_ns, 1).transpose();
			dmns(i + 1, j).boundary_north
			= gamma.block(start_i, 0, interface_size_ns, 1).transpose();
		}
	}

	// solve
	for (int i = 0; i < dmns.rows(); i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			dmns(i, j).solve();
		}
	}

	VectorXd diff(gamma.size());
	// get the difference on east-west interfaces
	for (int j = 0; j < dmns.cols() - 1; j++) {
		for (int i = 0; i < dmns.rows(); i++) {
			int start_i            = (j * dmns.rows() + i) * interface_size_ew;
			int left_dmns_last_col = dmns(i, j).grid.cols() - 1;
			diff.block(start_i, 0, interface_size_ew, 1)
			= dmns(i, j).u.col(left_dmns_last_col) + dmns(i, j + 1).u.col(0)
			  - 2 * gamma.block(start_i, 0, interface_size_ew, 1);
		}
	}
	// get the difference on north-south interfaces
	for (int i = 0; i < dmns.rows() - 1; i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			int start_i           = (i * dmns.cols() + j) * interface_size_ns + ns_start_i;
			int top_dmns_last_row = dmns(i, j).grid.rows() - 1;
			diff.block(start_i, 0, interface_size_ns, 1)
			= dmns(i, j).u.row(top_dmns_last_row).transpose() + dmns(i + 1, j).u.row(0).transpose()
			  - 2 * gamma.block(start_i, 0, interface_size_ns, 1);
		}
	}
	return diff;
}

SparseMatrix<double> generateMatrix(int m_x, int m_y, double h_x, double h_y)
{
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
	args::ArgumentParser  parser("");
	args::HelpFlag        help(parser, "help", "Display this help menu", {'h', "help"});
	args::Positional<int> d_x(parser, "d_x", "number of domains in the x direction");
	args::Positional<int> d_y(parser, "d_y", "number of domains in the y direction");
	args::Positional<int> n_x(parser, "n_x", "number of cells in the x direction, in each domain");
	args::Positional<int> n_y(parser, "n_y", "number of cells in the y direction, in each domain");
	args::ValueFlag<string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                            {'m'});
	args::ValueFlag<string> f_s(parser, "solution filename", "the file to write the solution to",
	                            {'s'});
	args::Flag f_cg(parser, "cg", "use conjugate gradient for solving gamma values", {"cg"});

	if (argc < 5) {
		std::cout << parser;
		return 0;
	}
	try {
		parser.ParseCLI(argc, argv);
	} catch (args::Help) {
		std::cout << parser;
		return 0;
	} catch (args::ParseError e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	} catch (args::ValidationError e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	int    num_domains_x  = args::get(d_x);
	int    num_domains_y  = args::get(d_y);
	int    domain_width_x = args::get(n_x);
	int    domain_width_y = args::get(n_y);
	int    m_x            = domain_width_x * num_domains_x;
	int    m_y            = domain_width_y * num_domains_y;
	double x_start        = 0.0;
	double x_end          = 1.0;
	double h_x            = (x_end - x_start) / m_x;
	double h_y            = (x_end - x_start) / m_y;

	string save_solution_file = "";
	if (f_s) {
		save_solution_file = args::get(f_s);
	}

	string save_matrix_file = "";
	if (f_m) {
		save_matrix_file = args::get(f_m);
	}

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
	DomainMatrix dmns(num_domains_y, num_domains_x);
	for (int i = 0; i < num_domains_y; i++) {
		for (int j = 0; j < num_domains_x; j++) {
			SparseMatrix<double> A       = generateMatrix(domain_width_x, domain_width_y, h_x, h_y);
			int                  start_i = i * domain_width_y;
			int                  start_j = j * domain_width_x;
			MatrixXd sub_grid = G.block(start_i, start_j, domain_width_y, domain_width_x);
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
		dmns(i, 0).boundary_west                 = boundary_west.block(start_i, 0, length, 1);
		dmns(i, num_domains_x - 1).boundary_east = boundary_east.block(start_i, 0, length, 1);
	}

	// create gamma array
	VectorXd gamma = VectorXd::Zero(m_y * (num_domains_x - 1) + m_x * (num_domains_y - 1));
	double   condition_number = 0.0;
	if (num_domains_x > 1 || num_domains_y > 1) {
		// get the b vector
		VectorXd b = solveWithInterfaceValues(dmns, gamma);

		if (f_cg) {
			FunctionWrapper F(&dmns, gamma.size(), b);
			Eigen::ConjugateGradient<FunctionWrapper, Eigen::Lower | Eigen::Upper,
			                         Eigen::IdentityPreconditioner>
			cg(F);
            cg.setTolerance(10e-10);
			// cg.compute(F);
			gamma = -1 * cg.solve(b);
			std::cout << "CG: Number of iterations: " << cg.iterations() << "\n";
		} else {
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

			if (save_matrix_file != "") {
				// print out solution
				ofstream out_file(save_matrix_file);
				out_file.precision(20);
				out_file << scientific;
				out_file << A << "\n";
				out_file.close();
			}
		}
	}

	// do one last solve
	solveWithInterfaceValues(dmns, gamma);

	// form the complete grid
	MatrixXd SOL(m_y, m_x);
	for (int i = 0; i < num_domains_y; i++) {
		for (int j = 0; j < num_domains_x; j++) {
			int start_i = i * domain_width_y;
			int start_j = j * domain_width_x;
			SOL.block(start_i, start_j, domain_width_y, domain_width_x) = dmns(i, j).u;
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

