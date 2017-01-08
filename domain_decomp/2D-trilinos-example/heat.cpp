#include "FunctionWrapper.h"
#include "MyTypeDefs.h"
#include "args.h"
#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>
#include <Epetra_CombineMode.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Vector.h>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
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

int main(int argc, char *argv[])
{
	//    using Belos::FuncWrap;
	using Teuchos::RCP;
	using Teuchos::rcp;

	MPI_Init(&argc, &argv);
	Epetra_MpiComm Comm(MPI_COMM_WORLD);

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);
	// comunicator for subdomain
	// for now, subdomains are going to be single threaded
	MPI_Comm subdomain_comm_raw;
	MPI_Comm_split(MPI_COMM_WORLD, my_global_rank, 0, &subdomain_comm_raw);
	Epetra_MpiComm subdomain_comm(subdomain_comm_raw);

	// parse input
	args::ArgumentParser  parser("");
	args::HelpFlag        help(parser, "help", "Display this help menu", {'h', "help"});
	args::Positional<int> d_x(parser, "d_x", "number of domains in the x direction");
	args::Positional<int> d_y(parser, "d_y", "number of domains in the y direction");
	args::Positional<int> n_x(parser, "n_x", "number of cells in the x direction, in each domain");
	args::Positional<int> n_y(parser, "n_y", "number of cells in the y direction, in each domain");

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
	int    domain_x      = my_global_rank % num_domains_x;
	int    domain_y      = my_global_rank / num_domains_x;

	if (num_domains_x * num_domains_y != num_procs) {
		std::cerr << "number of domains must be equal to the number of processes\n";
		return 1;
	}

	// create a map and matrix
	RCP<map_type> Map = rcp(new map_type(nx * ny, 0, subdomain_comm));

	// Generate RHS vector
	RCP<vector_type> f     = rcp(new vector_type(*Map, 1, false));
	RCP<vector_type> exact = rcp(new vector_type(*Map, 1, false));

	// Use local indices to access the entries of f_data.
	for (int yi = 0; yi < ny; yi++) {
		for (int xi = 0; xi < nx; xi++) {
			int    index_x            = domain_x * nx + xi;
			int    index_y            = domain_y * ny + yi;
			double x                  = h_x / 2.0 + 1.0 * index_x / (nx * num_domains_x);
			double y                  = h_y / 2.0 + 1.0 * index_y / (ny * num_domains_y);
			(*f)[0][yi * nx + xi]     = ffun(x, y);
			(*exact)[0][yi * nx + xi] = gfun(x, y);
			if (index_x == 0) {
				(*f)[0][yi * nx + xi] += -2.0 / (h_x * h_x) * gfun(0.0, y);
			}
			if (index_x == num_domains_x * nx - 1) {
				(*f)[0][yi * nx + xi] += -2.0 / (h_x * h_x) * gfun(1.0, y);
			}
			if (index_y == 0) {
				(*f)[0][yi * nx + xi] += -2.0 / (h_y * h_y) * gfun(x, 0.0);
			}
			if (index_y == num_domains_y * ny - 1) {
				(*f)[0][yi * nx + xi] += -2.0 / (h_y * h_y) * gfun(x, 1.0);
			}
		}
	}

	// create a domain
	Domain d(f, nx, ny, h_x, h_y);

	// Detertime on which sides this domain has a neighboring domain
	// and calculate how many interface points this domain is going to have
	int num_interface_points = 0;
	// north
	if (domain_y != num_domains_y - 1) {
		num_interface_points += nx;
		d.has_north = true;
	}
	// east
	if (domain_x != num_domains_x - 1) {
		num_interface_points += ny;
		d.has_east = true;
	}
	// south
	if (domain_y != 0) {
		num_interface_points += nx;
		d.has_south = true;
	}
	// west
	if (domain_x != 0) {
		num_interface_points += ny;
		d.has_west = true;
	}

	// Now Calculate the global indicies for thos interface points
	std::vector<int> global_i(num_interface_points);
	RCP<vector_type> u          = d.u;
	int              ns_start_i = num_domains_y * (num_domains_x - 1) * ny;
	int              curr_i     = 0;
	if (d.has_north) {
		int curr_global_i = (domain_y * num_domains_x + domain_x) * nx + ns_start_i;
		for (int i = 0; i < nx; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	if (d.has_east) {
		int curr_global_i = (domain_x * num_domains_y + domain_y) * ny;
		for (int i = 0; i < ny; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	if (d.has_south) {
		int curr_global_i = ((domain_y - 1) * num_domains_x + domain_x) * nx + ns_start_i;
		for (int i = 0; i < nx; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	if (d.has_west) {
		int curr_global_i = ((domain_x - 1) * num_domains_y + domain_y) * ny;
		for (int i = 0; i < ny; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}

	// Now that the global indices have been calculated, we can create a map for the interface
	// points
	int num_global_elements
	= nx * num_domains_x * (num_domains_y - 1) + ny * num_domains_y * (num_domains_x - 1);
	if (num_domains_x * num_domains_y == 1) {
		// this is a special case for when there is only one domain
		d.domain_map = rcp(new map_type(1, 0, Comm));
	} else {
		d.domain_map = rcp(new map_type(-1, num_interface_points, &global_i[0], 0, Comm));
	}

	// Create a map that will be used in the iterative solver
	RCP<map_type> diff_map = rcp(new map_type(num_global_elements, 0, Comm));

	// Create the gamma and diff vectors
	RCP<vector_type> gamma = rcp(new vector_type(*diff_map, 1));
	RCP<vector_type> diff  = rcp(new vector_type(*diff_map, 1));

	if (num_domains_x * num_domains_y != 1) {
		// do iterative solve

		// Get the b vector
		RCP<vector_type> b = rcp(new vector_type(*diff_map, 1));
		d.solveWithInterface(*gamma, *b);

		// Create a function wrapper
		RCP<FuncWrap> wrapper = rcp(new FuncWrap(b, &d));

		// Create linear problem for the Belos solver
		Belos::LinearProblem<double, vector_type, FuncWrap> problem(wrapper, gamma, b);
		problem.setProblem();

		// Set the parameters
		Teuchos::ParameterList belosList;
		belosList.set("Block Size", 1);
		belosList.set("Maximum Iterations", 1000);
		belosList.set("Convergence Tolerance", 10e-10);
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
	d.solveWithInterface(*gamma, *diff);

	// Calcuate error
	map_type      err_map(-1, 1, 0, Comm);
	Epetra_Vector exact_norm(err_map);
	Epetra_Vector diff_norm(err_map);
	exact->Norm2(&exact_norm[0]);
	{
		// Use local indices to access the entries of f_data.
		const int localLength      = f->MyLength();
		int       NumMyElements    = Map->NumMyElements();
		int *     MyGlobalElements = 0;
		Map->MyGlobalElementsPtr(MyGlobalElements);
		for (int i = 0; i < NumMyElements; i++) {
			(*exact)[0][i] -= (*u)[0][i];
		}
	}
	exact->Norm2(&diff_norm[0]);
	double global_diff_norm;
	double global_exact_norm;
	diff_norm.Norm2(&global_diff_norm);
	exact_norm.Norm2(&global_exact_norm);
	if (my_global_rank == 0) {
        std::cout << std::scientific;
        std::cout.precision(13);
		std::cout << "Error: " << global_diff_norm / global_exact_norm << "\n";
	}

	MPI_Finalize();
	return 0;
}
