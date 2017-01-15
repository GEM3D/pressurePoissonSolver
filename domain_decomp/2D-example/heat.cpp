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
	DomainMatrix *     dmns;
	vector<Interface> *interfaces;
	VectorXd           b;
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
	FunctionWrapper(DomainMatrix *dmns, vector<Interface> *interfaces, int size, VectorXd b)
	{
		this->dmns       = dmns;
		this->interfaces = interfaces;
		this->b          = b;
		my_rows          = size;
	}

	void attachMyMatrix(const SparseMatrix<double> &mat) { mp_mat = &mat; }
	const SparseMatrix<double>                      my_matrix() const { return *mp_mat; }
};
VectorXd solveWithInterfaceValues(DomainMatrix &dmns, vector<Interface> &interfaces,
                                  VectorXd &gamma);
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
		VectorXd result  = lhs.b - solveWithInterfaceValues(*lhs.dmns, *lhs.interfaces, rhscopy);
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
VectorXd solveWithInterfaceValues(DomainMatrix &dmns, vector<Interface> &interfaces,
                                  VectorXd &gamma)
{
	for (Interface &i : interfaces) {
		if (i.dir == Interface::axis::y) {
			i.gamma = gamma.block(i.start_index, 0, i.size, 1);
		} else {
			i.gamma = gamma.block(i.start_index, 0, i.size, 1).transpose();
		}
	}

	// solve
	for (int i = 0; i < dmns.rows(); i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			dmns(i, j).solve();
		}
	}

	VectorXd diff(gamma.size());
	for (Interface &i : interfaces) {
		diff.block(i.start_index, 0, i.size, 1) = i.getDiff();
	}
	return diff;
}

vector<Interface> createAndLinkInterfaces(DomainMatrix &dmns)
{
	int               num_interfaces_ns = dmns.cols() * (dmns.rows() - 1);
	int               num_interfaces_ew = (dmns.cols() - 1) * dmns.rows();
	int               num_interfaces    = num_interfaces_ns + num_interfaces_ew;
	int               ns_interface_size = dmns(0, 0).u.cols();
	int               ew_interface_size = dmns(0, 0).u.rows();
	vector<Interface> interfaces(num_interfaces);
	// fill the array of interfaces
	int curr_i = 0;
	for (int j = 0; j < dmns.cols() - 1; j++) {
		for (int i = 0; i < dmns.rows(); i++) {
			int start_i = (j * dmns.rows() + i) * ew_interface_size;
			interfaces[curr_i]        = Interface(start_i, ew_interface_size, Interface::axis::y);
			Interface &curr_interface = interfaces[curr_i];
			// link to domain on left
			curr_interface.left = &dmns(i, j);
			dmns(i, j).east = &curr_interface;
			// link to domain on right
			curr_interface.right = &dmns(i, j + 1);
			dmns(i, j+1).west = &curr_interface;
			curr_i++;
		}
	}

	int ns_start_i = dmns.rows() * (dmns.cols() - 1) * ew_interface_size;
	for (int i = 0; i < dmns.rows() - 1; i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			int start_i = (i * dmns.cols() + j) * ns_interface_size + ns_start_i;
			interfaces[curr_i]        = Interface(start_i, ns_interface_size, Interface::axis::x);
			Interface &curr_interface = interfaces[curr_i];
			// link to domain on left
			curr_interface.left = &dmns(i, j);
			dmns(i, j).north = &curr_interface;
			// link to domain on right
			curr_interface.right = &dmns(i + 1, j);
			dmns(i+1, j).south = &curr_interface;
			curr_i++;
		}
	}
	return interfaces;
}
void graphAssistedMatrixFormation(DomainMatrix &dmns, MatrixXd &A, VectorXd &b) {}
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
	args::Flag f_gp(parser, "graph", "use a graph when forming the matrix", {"graph"});

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
		boundary_east(i) = gfun(x_end, y(i));
		boundary_west(i) = gfun(x_start, y(i));
	}

	RowVectorXd boundary_north(m_x);
	RowVectorXd boundary_south(m_x);
	for (int j = 0; j < m_x; j++) {
		boundary_north(j) = gfun(x(j), x_end);
		boundary_south(j) = gfun(x(j), x_start);
	}

	// set outer boundaries
	// north
	G.row(G.rows() - 1) -= 2 / ((h_y) * (h_y)) * boundary_north;
	// east
	G.col(G.cols() - 1) -= 2 / ((h_x) * (h_x)) * boundary_east;
	// south
	G.row(0) -= 2 / ((h_y) * (h_y)) * boundary_south;
	// west
	G.col(0) -= 2 / ((h_x) * (h_x)) * boundary_west;

	// Generate Domains
	DomainMatrix dmns(num_domains_y, num_domains_x);
	for (int i = 0; i < num_domains_y; i++) {
		for (int j = 0; j < num_domains_x; j++) {
			int      start_i  = i * domain_width_y;
			int      start_j  = j * domain_width_x;
			MatrixXd sub_grid = G.block(start_i, start_j, domain_width_y, domain_width_x);
			dmns(i, j) = Domain(sub_grid, h_x, h_y);
		}
	}

	// Generate Interfaces
	vector<Interface> interfaces = createAndLinkInterfaces(dmns);
	// create gamma array
	VectorXd gamma = VectorXd::Zero(m_y * (num_domains_x - 1) + m_x * (num_domains_y - 1));
	double   condition_number = 0.0;
	if (num_domains_x > 1 || num_domains_y > 1) {
		// get the b vector
		VectorXd b = solveWithInterfaceValues(dmns, interfaces, gamma);

		if (f_cg) {
			FunctionWrapper F(&dmns, &interfaces, gamma.size(), b);
			Eigen::ConjugateGradient<FunctionWrapper, Eigen::Lower | Eigen::Upper,
			                         Eigen::IdentityPreconditioner>
			cg(F);
			cg.setTolerance(1e-10);
			// cg.compute(F);
			gamma = cg.solve(b);
			std::cout << "CG: Number of iterations: " << cg.iterations() << "\n";
		} else {
			// create a matrix
			MatrixXd A(gamma.size(), gamma.size());
			// get the columns of the matrix
			for (int i = 0; i < gamma.size(); i++) {
				gamma(i) = 1;
				A.col(i) = b - solveWithInterfaceValues(dmns, interfaces, gamma);
				gamma(i) = 0;
			}

			cout << "Number of solves to form " << gamma.size() << "x" << gamma.size()
			     << " Matrix: " << gamma.size() << "\n";

			// solve for gamme values
			FullPivLU<MatrixXd> lu(A);
			gamma            = lu.solve(b);
			condition_number = 1.0 / lu.rcond();

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
	solveWithInterfaceValues(dmns, interfaces, gamma);

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

