#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/SparseExtra>
#include <utility>

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

/**
 * @brief Create interfaces, and set them up so that they are linked in a pre-defined order
 *
 * @param dmns the domains to use
 *
 * @return a vector of interface objects
 */
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
			interfaces[curr_i]        = Interface(curr_i,start_i, ew_interface_size, Interface::axis::y);
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
			interfaces[curr_i]        = Interface(curr_i,start_i, ns_interface_size, Interface::axis::x);
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
/**
 * @brief Create interfaces, and set them up using a Breadth-First Search
 *
 * @param dmns the domains to use
 *
 * @return a vector of interface objects
 */
vector<Interface> createAndLinkInterfacesBFS(DomainMatrix &dmns)
{
	int               num_interfaces_ns = dmns.cols() * (dmns.rows() - 1);
	int               num_interfaces_ew = (dmns.cols() - 1) * dmns.rows();
	int               num_interfaces    = num_interfaces_ns + num_interfaces_ew;
	vector<Interface> interfaces(num_interfaces);
	// fill the array of interfaces
	deque<Domain *> queue;
	queue.push_back(&dmns(0, 0));
	set<Domain *> visited;
	int           curr_i       = 0;
	int           curr_gamma_i = 0;
	while (!queue.empty()) {
		Domain &curr_domain = *queue.front();
		queue.pop_front();
		visited.insert(&curr_domain);
		if (curr_domain.north_domain != nullptr && visited.count(curr_domain.north_domain) == 0) {
			Domain &nbr            = *curr_domain.north_domain;
			int     interface_size = curr_domain.u.cols();
			interfaces[curr_i]
			= Interface(curr_i, curr_gamma_i, interface_size, Interface::axis::x);
            interfaces[curr_i].left = &curr_domain;
            curr_domain.north = &interfaces[curr_i];
            interfaces[curr_i].right = &nbr;
            nbr.south = &interfaces[curr_i];
            curr_i++;
            curr_gamma_i+=interface_size;
			if (std::find(queue.begin(), queue.end(), &nbr) == queue.end()) {
				queue.push_back(&nbr);
			}
		}
		if (curr_domain.east_domain != nullptr && visited.count(curr_domain.east_domain) == 0) {
			Domain &nbr            = *curr_domain.east_domain;
			int     interface_size = curr_domain.u.rows();
			interfaces[curr_i]
			= Interface(curr_i, curr_gamma_i, interface_size, Interface::axis::y);
            interfaces[curr_i].left = &curr_domain;
            curr_domain.east = &interfaces[curr_i];
            interfaces[curr_i].right = &nbr;
            nbr.west = &interfaces[curr_i];
            curr_i++;
            curr_gamma_i+=interface_size;
			if (std::find(queue.begin(), queue.end(), &nbr) == queue.end()) {
				queue.push_back(&nbr);
			}
		}
		if (curr_domain.south_domain != nullptr && visited.count(curr_domain.south_domain) == 0) {
			Domain &nbr            = *curr_domain.south_domain;
			int     interface_size = curr_domain.u.cols();
			interfaces[curr_i]
			= Interface(curr_i, curr_gamma_i, interface_size, Interface::axis::x);
            interfaces[curr_i].left = &nbr;
            nbr.north = &interfaces[curr_i];
            interfaces[curr_i].right = &curr_domain;
            curr_domain.south = &interfaces[curr_i];
            curr_i++;
            curr_gamma_i+=interface_size;
			if (std::find(queue.begin(), queue.end(), &nbr) == queue.end()) {
				queue.push_back(&nbr);
			}
		}
		if (curr_domain.west_domain != nullptr && visited.count(curr_domain.west_domain) == 0) {
			Domain &nbr            = *curr_domain.west_domain;
			int     interface_size = curr_domain.u.rows();
			interfaces[curr_i]
			= Interface(curr_i, curr_gamma_i, interface_size, Interface::axis::y);
            interfaces[curr_i].left = &nbr;
            nbr.east = &interfaces[curr_i];
            interfaces[curr_i].right = &curr_domain;
            curr_domain.west = &interfaces[curr_i];
            curr_i++;
            curr_gamma_i+=interface_size;
			if (std::find(queue.begin(), queue.end(), &nbr) == queue.end()) {
				queue.push_back(&nbr);
			}
		}
	}
	return interfaces;
}
void graphAssistedMatrixFormation(DomainMatrix &dmns, vector<Interface> &interfaces, SparseMatrix<double> &A,
                                  VectorXd &b)
{
	using namespace boost;
	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
	// typedef std::pair<int, int> Edge;
	// typedef graph_traits<Graph>::vertex_descriptor  vertex_descriptor;
	typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
	typedef property_map<Graph, vertex_index_t>::const_type vertex_index_map;

	// create some enums
	enum { NORTH, EAST, SOUTH, WEST };
	enum { LEFT, RIGHT };

	const int num_vertices = interfaces.size();
	Graph     graph(num_vertices);

	// populate graph with edges
	for (Interface &curr_iface : interfaces) {
		const int u = curr_iface.my_id;
		for (int i = LEFT; i <= RIGHT; i++) {
			// get domain
			Domain *domain_ptr;
			if (i == LEFT) {
				domain_ptr = curr_iface.left;
			} else {
				domain_ptr = curr_iface.right;
			}
			Domain &domain = *domain_ptr;

			for (int i = NORTH; i <= WEST; i++) {
				// get neighbor interface
				Interface *nbr_iface;
				switch (i) {
					case NORTH:
						nbr_iface = domain.north;
						break;
					case EAST:
						nbr_iface = domain.east;
						break;
					case SOUTH:
						nbr_iface = domain.south;
						break;
					case WEST:
						nbr_iface = domain.west;
						break;
				}

				// add an edge for that interface
				if (nbr_iface != nullptr && nbr_iface != &curr_iface) {
					int v = nbr_iface->my_id;
					add_edge(u, v, graph);
				}
			}
		}
	}

	std::vector<size_t> color_vec(num_vertices);
	iterator_property_map<size_t *, vertex_index_map> color(&color_vec.front(),
	                                                        get(vertex_index, graph));
	vertices_size_type num_colors = sequential_vertex_coloring(graph, color);
	cout << "Number of colors: " << num_colors << "\n";
	cout << "Colors: \n";
	for (size_t c : color_vec) {
		cout << c << " ";
	}
	cout << "\n";
	// get b vector
	VectorXd gamma = VectorXd::Zero(b.size());
	b              = solveWithInterfaceValues(dmns, interfaces, gamma);
	// save the edge results for zero, and get max interface size
	int max_size = 0;
	for (Interface &i : interfaces) {
		i.setZero();
		if (i.size > max_size) {
			max_size = i.size;
		}
	}

	typedef Eigen::Triplet<double> Trip;
	std::deque<Trip>               tripletList;
	// start forming matrix
	int num_solves = 0;
	for (size_t c = 0; c < num_colors; c++) {
		// get list of vertices with color c
		list<Interface *> vertices;
		for (size_t i = 0; i < color_vec.size(); i++) {
			if (color_vec[i] == c) {
				vertices.push_back(&interfaces[i]);
			}
		}
		// set each value to one and solve
		for (int index = 0; index < max_size; index++) {
			num_solves++;
			// set value to one
			for (Interface *curr_iface : vertices) {
				if (index < curr_iface->size) {
					curr_iface->gamma(index) = 1;
				}
			}

			// solve
			for (int i = 0; i < dmns.rows(); i++) {
				for (int j = 0; j < dmns.cols(); j++) {
					dmns(i, j).solve();
				}
			}

			// fill matrix values
			for (Interface *curr_iface : vertices) {
				if (index < curr_iface->size) {
					// set diagonal
					{
						int      iface_i = curr_iface->start_index;
						int      j       = iface_i + index;
						VectorXd diff    = b.block(curr_iface->start_index, 0, curr_iface->size, 1)
						                - (VectorXd) curr_iface->getDiff();
						for (int curr_i = 0; curr_i < curr_iface->size; curr_i++) {
							int i = curr_i + iface_i;
							tripletList.push_back(Trip(i, j, diff(curr_i)));
						}
					}

					// set off-diagonal entries
					for (int i = LEFT; i <= RIGHT; i++) {
						// get domain
						Domain *domain_ptr;
						if (i == LEFT) {
							domain_ptr = curr_iface->left;
						} else {
							domain_ptr = curr_iface->right;
						}
						Domain &domain = *domain_ptr;

						for (int i = NORTH; i <= WEST; i++) {
							// get neighbor interface
							Interface *nbr_iface;
							switch (i) {
								case NORTH:
									nbr_iface = domain.north;
									break;
								case EAST:
									nbr_iface = domain.east;
									break;
								case SOUTH:
									nbr_iface = domain.south;
									break;
								case WEST:
									nbr_iface = domain.west;
									break;
							}

							// add entries for that interface
							if (nbr_iface != nullptr && nbr_iface != curr_iface) {
								// get the difference
								VectorXd diff
								= b.block(nbr_iface->start_index, 0, nbr_iface->size, 1);
								if (i == NORTH || i == EAST) {
									diff -= (VectorXd) nbr_iface->getZeroDiffFromLeft();
								} else {
									diff -= (VectorXd) nbr_iface->getZeroDiffFromRight();
								}

								// fill matrix values
								int iface_i = nbr_iface->start_index;
								int j       = curr_iface->start_index + index;
								for (int curr_i = 0; curr_i < nbr_iface->size; curr_i++) {
									int i = curr_i + iface_i;
									tripletList.push_back(Trip(i, j, diff(curr_i)));
								}
							}
						}
					}
				}
			}

			// set value back to zero
			for (Interface *curr_iface : vertices) {
				if (index < curr_iface->size) {
					curr_iface->gamma(index) = 0;
				}
			}
		}
	}
	// form a from the list of triplets
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	cout << "Number of solves to form " << A.rows() << "x" << A.cols() << " Matrix: " << num_solves
	     << "\n";
}

void linkDomains(DomainMatrix &dmns)
{
	for (int i = 0; i < dmns.rows(); i++) {
		for (int j = 0; j < dmns.cols(); j++) {
			if (i != dmns.rows() - 1) {
				dmns(i, j).north_domain = &dmns(i + 1, j);
			}
			if (j != dmns.cols() - 1) {
				dmns(i, j).east_domain = &dmns(i, j + 1);
			}
			if (i != 0) {
				dmns(i, j).south_domain = &dmns(i - 1, j);
			}
			if (j != 0) {
				dmns(i, j).west_domain = &dmns(i, j - 1);
			}
		}
	}
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
	args::Flag f_gp(parser, "graph", "use a graph when forming the matrix", {"graph"});
	args::Flag f_bfs(
	parser, "bfs",
	"use a breadth-first search when determining mapping of interfaces to gamma vectors", {"bfs"});

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
	linkDomains(dmns);

	// Generate Interfaces
	vector<Interface> interfaces;
	if (f_bfs) {
		interfaces = createAndLinkInterfacesBFS(dmns);
	} else {
		interfaces = createAndLinkInterfaces(dmns);
	}

	// create gamma array
	VectorXd gamma = VectorXd::Zero(m_y * (num_domains_x - 1) + m_x * (num_domains_y - 1));
	double   condition_number = 0.0;
	if (num_domains_x > 1 || num_domains_y > 1) {
		if (f_cg) {
			// Solve with a function wrapper
			// get the b vector
			VectorXd        b = solveWithInterfaceValues(dmns, interfaces, gamma);
			FunctionWrapper F(&dmns, &interfaces, gamma.size(), b);
			Eigen::ConjugateGradient<FunctionWrapper, Eigen::Lower | Eigen::Upper,
			                         Eigen::IdentityPreconditioner>
			cg(F);
			cg.setTolerance(1e-10);
			// cg.compute(F);
			gamma = cg.solve(b);
			std::cout << "CG: Number of iterations: " << cg.iterations() << "\n";
		} else if (f_gp) {
			// quickly form the matrix and then use CG to solve
			SparseMatrix<double> A(gamma.size(), gamma.size());
			VectorXd             b(gamma.size());

			// form the matrix
			graphAssistedMatrixFormation(dmns, interfaces, A, b);

			// solve
			Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
			                         Eigen::IdentityPreconditioner>
			cg(A);
			cg.setTolerance(1e-10);
			gamma = cg.solve(b);
			std::cout << "CG: Number of iterations: " << cg.iterations() << "\n";
			if (save_matrix_file != "") {
				saveMarket(A, save_matrix_file);
			}
		} else {
			// slowly form the matrix and then solve with LU decomposition
			// get the b vector
			VectorXd b = solveWithInterfaceValues(dmns, interfaces, gamma);
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

			// solve for gamma values
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

