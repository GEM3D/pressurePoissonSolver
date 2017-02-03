#include "FunctionWrapper.h"
#include "MyTypeDefs.h"
#include "args.h"
#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <unistd.h>

using Teuchos::RCP;
using Teuchos::rcp;

// the functions that we are using
double ffun(double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); }
double gfun(double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); }
// =========== //
// main driver //
// =========== //

using namespace std;
DomainCollection::DomainCollection(int low, int high, int nx, int ny, int d_x, int d_y, double h_x,
                                   double h_y, RCP<const Teuchos::Comm<int>> comm)
{
	// cerr<< "Low:  " << low << "\n";
	// cerr<< "High: " << high << "\n";
	this->comm = comm;
	this->nx   = nx;
	this->ny   = ny;
	this->h_x  = h_x;
	this->h_y  = h_y;
	domains    = map<int, Domain *>();
	for (int i = low; i <= high; i++) {
		int domain_x = i % d_y;
		int domain_y = i / d_x;

		// Generate RHS vector
		std::valarray<double> f(nx * ny);
		std::valarray<double> exact(nx * ny);
		// Use local indices to access the entries of f_data.
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				int    index_x      = domain_x * nx + xi;
				int    index_y      = domain_y * ny + yi;
				double x            = h_x / 2.0 + 1.0 * index_x / (nx * d_x);
				double y            = h_y / 2.0 + 1.0 * index_y / (ny * d_y);
				f[yi * nx + xi]     = ffun(x, y);
				exact[yi * nx + xi] = gfun(x, y);
				if (index_x == 0) {
					f[yi * nx + xi] += -2.0 / (h_x * h_x) * gfun(0.0, y);
				}
				if (index_x == d_x * nx - 1) {
					f[yi * nx + xi] += -2.0 / (h_x * h_x) * gfun(1.0, y);
				}
				if (index_y == 0) {
					f[yi * nx + xi] += -2.0 / (h_y * h_y) * gfun(x, 0.0);
				}
				if (index_y == d_y * ny - 1) {
					f[yi * nx + xi] += -2.0 / (h_y * h_y) * gfun(x, 1.0);
				}
			}
		}
		// create a domain
		Domain *d_ptr = new Domain(f, exact, nx, ny, h_x, h_y);
		Domain &d     = *d_ptr;

		// determine its neighbors
		// north
		int ns_start_i = d_y * (d_x - 1) * ny;
		if (domain_y != d_y - 1) {
			int nbr_y        = domain_y + 1;
			d.nbr_north      = nbr_y * d_x + domain_x;
			d.global_i_north = (domain_y * d_x + domain_x) * nx + ns_start_i;
		}
		// east
		if (domain_x != d_x - 1) {
			int nbr_x       = domain_x + 1;
			d.nbr_east      = domain_y * d_x + nbr_x;
			d.global_i_east = (domain_x * d_x + domain_y) * ny;
		}
		// south
		if (domain_y != 0) {
			int nbr_y        = domain_y - 1;
			d.nbr_south      = nbr_y * d_x + domain_x;
			d.global_i_south = ((domain_y - 1) * d_x + domain_x) * nx + ns_start_i;
		}
		// west
		if (domain_x != 0) {
			int nbr_x       = domain_x - 1;
			d.nbr_west      = domain_y * d_x + nbr_x;
			d.global_i_west = ((domain_x - 1) * d_x + domain_y) * ny;
		}
		domains[i] = d_ptr;
	}

	// create map for domains
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	queue.push_back(low);
	enqueued.insert(low);
	int         curr_i = 0;
	vector<int> global;
	global.reserve((high - low + 1) * (2 * nx + 2 * ny));
	while (!queue.empty()) {
		int     curr = queue.front();
		Domain &d    = *domains[curr];
		queue.pop_front();
		visited.insert(curr);
		if (d.nbr_north != -1 && visited.count(d.nbr_north) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr       = *domains.at(d.nbr_north);
				nbr.local_i_south = curr_i;
				if (enqueued.count(d.nbr_north) == 0) {
					queue.push_back(d.nbr_north);
					enqueued.insert(d.nbr_north);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_north = curr_i;
			for (int i = 0; i < nx; i++) {
				global.push_back(d.global_i_north + i);
			}
			curr_i += nx;
		}
		if (d.nbr_east != -1 && visited.count(d.nbr_east) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr      = *domains.at(d.nbr_east);
				nbr.local_i_west = curr_i;
				if (enqueued.count(d.nbr_east) == 0) {
					queue.push_back(d.nbr_east);
					enqueued.insert(d.nbr_east);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_east = curr_i;
			for (int i = 0; i < ny; i++) {
				global.push_back(d.global_i_east + i);
			}
			curr_i += ny;
		}
		if (d.nbr_south != -1 && visited.count(d.nbr_south) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr       = *domains.at(d.nbr_south);
				nbr.local_i_north = curr_i;
				if (enqueued.count(d.nbr_south) == 0) {
					queue.push_back(d.nbr_south);
					enqueued.insert(d.nbr_south);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_south = curr_i;
			for (int i = 0; i < nx; i++) {
				global.push_back(d.global_i_south + i);
			}
			curr_i += nx;
		}
		if (d.nbr_west != -1 && visited.count(d.nbr_west) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr      = *domains.at(d.nbr_west);
				nbr.local_i_east = curr_i;
				if (enqueued.count(d.nbr_west) == 0) {
					queue.push_back(d.nbr_west);
					enqueued.insert(d.nbr_west);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_west = curr_i;
			for (int i = 0; i < ny; i++) {
				global.push_back(d.global_i_west + i);
			}
			curr_i += ny;
		}
	}
	// Now that the global indices have been calculated, we can create a map for the interface
	// points
	if (d_x * d_y == 1) {
		// this is a special case for when there is only one domain
		collection_map = Teuchos::rcp(new map_type(1, 0, comm));
	} else {
		int num_global_elements = nx * d_x * (d_y - 1) + ny * d_y * (d_x - 1);
		collection_map
		= Teuchos::rcp(new map_type(num_global_elements, &global[0], curr_i, 0, this->comm));
	}
}
void DomainCollection::solveWithInterface(const vector_type &gamma, vector_type &diff)
{
	// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
	// if(has_east)std::cout << "Gamma begin\n";
	// gamma.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.update(1, gamma, 0);
	Tpetra::Import<> importer(diff.getMap(), collection_map);
	Tpetra::Export<> exporter(collection_map, diff.getMap());
	vector_type      local_gamma(collection_map, 1);
	vector_type      local_diff(collection_map, 1);
	local_gamma.doImport(diff, importer, Tpetra::CombineMode::INSERT);

	// solve over domains on this proc
	for (auto &p : domains) {
		p.second->solveWithInterface(local_gamma, local_diff);
	}

	// export diff vector

	// if(has_east)std::cout <<"LOCAL AFEr\n";
	// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.scale(0);
	// if(has_east)std::cout <<"DIFFBEFORE\n";
	// diff.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.doExport(local_diff, importer, Tpetra::CombineMode::ADD);
	// if(has_east)std::cout <<"DIFF\n";
	// diff.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	// if(has_east)std::cout <<"GAMMA\n";
	// gamma.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.update(-2, gamma, 1);
}
double DomainCollection::diffNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->diffNorm(), 2);
	}
	return sqrt(result);
}
double DomainCollection::exactNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->exactNorm(), 2);
	}
	return sqrt(result);
}
RCP<matrix_type> DomainCollection::formMatrix(RCP<map_type> map)
{
	// create domain for forming matrix
	std::valarray<double> f(nx * ny);
	Domain                d(f, f, nx, ny, h_x, h_y);
	d.nbr_north           = 1;
	d.nbr_east            = 1;
	d.nbr_south           = 1;
	d.nbr_west            = 1;
	d.boundary_north      = valarray<double>(nx);
	d.boundary_south      = valarray<double>(nx);
	d.boundary_east       = valarray<double>(ny);
	d.boundary_west       = valarray<double>(ny);
	int              size = max(nx, ny);
	RCP<matrix_type> A    = rcp(new matrix_type(map, size * 6));
	// north boundary
	for (int i = 0; i < nx; i++) {
		d.boundary_north[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_north != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_north + i;
				for (int i = 0; i < nx; i++) {
					row.push_back(d.u[nx * (ny - 1) + i]);
					global.push_back(d2.global_i_north + i);
				}
				row[i] -= 1;
				if (d2.nbr_east != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(d.u[(i + 1) * nx - 1]);
						global.push_back(d2.global_i_east + i);
					}
				}
				if (d2.nbr_south != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(d.u[i]);
						global.push_back(d2.global_i_south + i);
					}
				}
				if (d2.nbr_west != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(d.u[i * nx]);
						global.push_back(d2.global_i_west + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_north[i] = 0;
	}
	// east boundary
	for (int i = 0; i < ny; i++) {
		d.boundary_east[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_east != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_east + i;
				for (int i = 0; i < ny; i++) {
					row.push_back(d.u[(i + 1) * nx - 1]);
					global.push_back(d2.global_i_east + i);
				}
				row[i] -= 1;
				if (d2.nbr_north != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(d.u[nx * (ny - 1) + i]);
						global.push_back(d2.global_i_north + i);
					}
				}
				if (d2.nbr_south != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(d.u[i]);
						global.push_back(d2.global_i_south + i);
					}
				}
				if (d2.nbr_west != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(d.u[i * nx]);
						global.push_back(d2.global_i_west + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_east[i] = 0;
	}
	// south boundary
	for (int i = 0; i < nx; i++) {
		d.boundary_south[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_south != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_south + i;
				for (int i = 0; i < nx; i++) {
					row.push_back(d.u[i]);
					global.push_back(d2.global_i_south + i);
				}
				row[i] -= 1;
				if (d2.nbr_north != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(d.u[nx * (ny - 1) + i]);
						global.push_back(d2.global_i_north + i);
					}
				}
				if (d2.nbr_east != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(d.u[(i + 1) * nx - 1]);
						global.push_back(d2.global_i_east + i);
					}
				}
				if (d2.nbr_west != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(d.u[i * nx]);
						global.push_back(d2.global_i_west + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_south[i] = 0;
	}
	// west boundary
	for (int i = 0; i < ny; i++) {
		d.boundary_west[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_west != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_west + i;
				for (int i = 0; i < ny; i++) {
					row.push_back(d.u[i * nx]);
					global.push_back(d2.global_i_west + i);
				}
				row[i] -= 1;
				if (d2.nbr_north != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(d.u[nx * (ny - 1) + i]);
						global.push_back(d2.global_i_north + i);
					}
				}
				if (d2.nbr_east != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(d.u[(i + 1) * nx - 1]);
						global.push_back(d2.global_i_east + i);
					}
				}
				if (d2.nbr_south != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(d.u[i]);
						global.push_back(d2.global_i_south + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_west[i] = 0;
	}
	// transpose matrix and return
	A->fillComplete();
	return A;
}
int main(int argc, char *argv[])
{
	//    using Belos::FuncWrap;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using namespace std::chrono;

	MPI_Init(&argc, &argv);
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// parse input
	args::ArgumentParser  parser("");
	args::HelpFlag        help(parser, "help", "Display this help menu", {'h', "help"});
	args::Positional<int> d_x(parser, "d_x", "number of domains in the x direction");
	args::Positional<int> d_y(parser, "d_y", "number of domains in the y direction");
	args::Positional<int> n_x(parser, "n_x", "number of cells in the x direction, in each domain");
	args::Positional<int> n_y(parser, "n_y", "number of cells in the y direction, in each domain");
	args::ValueFlag<string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                            {'m'});

	if (argc < 5) {
		if (my_global_rank == 0) std::cout << parser;
		return 0;
	}
	try {
		parser.ParseCLI(argc, argv);
	} catch (args::Help) {
		if (my_global_rank == 0) std::cout << parser;
		return 0;
	} catch (args::ParseError e) {
		if (my_global_rank == 0) {
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
		}
		return 1;
	} catch (args::ValidationError e) {
		if (my_global_rank == 0) {
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
		}
		return 1;
	}
	// Set the number of discretization points in the x and y direction.
	int    num_domains_x = args::get(d_x);
	int    num_domains_y = args::get(d_y);
	int    nx            = args::get(n_x);
	int    ny            = args::get(n_y);
	double h_x           = 1.0 / (nx * num_domains_x);
	double h_y           = 1.0 / (ny * num_domains_y);

	if (num_domains_x * num_domains_y < num_procs) {
		std::cerr << "number of domains must be greater than or equal to the number of processes\n";
		return 1;
	}

	string save_matrix_file = "";
	if (f_m) {
		save_matrix_file = args::get(f_m);
	}

	int              total_domains = num_domains_x * num_domains_y;
	DomainCollection dc(total_domains * my_global_rank / num_procs,
	                    total_domains * (my_global_rank + 1) / num_procs - 1, nx, ny, num_domains_x,
	                    num_domains_y, h_x, h_y, comm);

	MPI_Barrier(MPI_COMM_WORLD);
	steady_clock::time_point iter_start = steady_clock::now();
	// Create a map that will be used in the iterative solver
	int num_global_elements
	= nx * num_domains_x * (num_domains_y - 1) + ny * num_domains_y * (num_domains_x - 1);
	RCP<map_type> diff_map = rcp(new map_type(num_global_elements, 0, comm));

	// Create the gamma and diff vectors
	RCP<vector_type> gamma = rcp(new vector_type(diff_map, 1));
	RCP<vector_type> diff  = rcp(new vector_type(diff_map, 1));

	if (num_domains_x * num_domains_y != 1) {
		// do iterative solve

		// Get the b vector
		RCP<vector_type> b = rcp(new vector_type(diff_map, 1));
		dc.solveWithInterface(*gamma, *b);

		// Create a function wrapper
		RCP<FuncWrap> wrapper = rcp(new FuncWrap(b, &dc));

		// Create linear problem for the Belos solver
		Belos::LinearProblem<double, vector_type, FuncWrap> problem(wrapper, gamma, b);
		problem.setProblem();

		// Set the parameters
		Teuchos::ParameterList belosList;
		belosList.set("Block Size", 1);
		belosList.set("Maximum Iterations", 1000);
		belosList.set("Convergence Tolerance", 1e-10);
		int verbosity = Belos::Errors + Belos::StatusTestDetails + Belos::Warnings
		                + Belos::TimingDetails + Belos::Debug;
		belosList.set("Verbosity", verbosity);
		Belos::OutputManager<double> my_om();

		// Create solver and solve
		RCP<Belos::SolverManager<double, vector_type, FuncWrap>> solver
		= rcp(new Belos::BlockCGSolMgr<double, vector_type, FuncWrap>(rcp(&problem, false),
		                                                              rcp(&belosList, false)));
		solver->solve();
	}

	// Do one last solve
	dc.solveWithInterface(*gamma, *diff);

	// Calcuate error
	RCP<map_type>    err_map = rcp(new map_type(-1, 1, 0, comm));
	Tpetra::Vector<> exact_norm(err_map);
	Tpetra::Vector<> diff_norm(err_map);

	exact_norm.getDataNonConst()[0] = dc.exactNorm();
	diff_norm.getDataNonConst()[0]  = dc.diffNorm();
	double global_diff_norm;
	double global_exact_norm;
	global_diff_norm  = diff_norm.norm2();
	global_exact_norm = exact_norm.norm2();
	MPI_Barrier(MPI_COMM_WORLD);
	duration<double> iter_time = steady_clock::now() - iter_start;
	if (my_global_rank == 0) {
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error: " << global_diff_norm / global_exact_norm << "\n";
		std::cout << std::defaultfloat;
		std::cout << "Time: " << iter_time.count() << "\n";
	}
	if (save_matrix_file != "") {
		MPI_Barrier(MPI_COMM_WORLD);
		steady_clock::time_point form_start = steady_clock::now();

		RCP<matrix_type> A = dc.formMatrix(diff_map);

		MPI_Barrier(MPI_COMM_WORLD);
		duration<double> form_time = steady_clock::now() - form_start;

		if (my_global_rank == 0) cout << "Matrix Formation Time: " << form_time.count() << "\n";

		MPI_Barrier(MPI_COMM_WORLD);
		steady_clock::time_point write_start = steady_clock::now();

		Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(save_matrix_file, A, "", "");

		MPI_Barrier(MPI_COMM_WORLD);
		duration<double> write_time = steady_clock::now() - write_start;
		if (my_global_rank == 0)
			cout << "Time to write matix to file: " << write_time.count() << "\n";
	}

	MPI_Finalize();
	return 0;
}
