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
#include <Teuchos_Comm.hpp>
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

// =========== //
// main driver //
// =========== //

using namespace std;
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
	args::ValueFlag<string> f_r(parser, "rhs filename", "the file to write the rhs vector to",
	                            {'r'});
	args::ValueFlag<string> f_p(parser, "preconditioner filename",
	                            "the file to write the preconditoiner to", {'p'});
	args::Flag f_wrapper(parser, "wrapper", "use a function wrapper", {"wrap"});
	args::Flag f_gauss(parser, "gauss", "solve gaussian function", {"gauss"});
	args::Flag f_prec(parser, "prec", "use block diagonal preconditioner", {"prec"});
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {'n', "neumann"});

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
	int    total_cells   = nx * num_domains_x * ny * num_domains_y;
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

	string save_rhs_file = "";
	if (f_r) {
		save_rhs_file = args::get(f_r);
	}

	string save_prec_file = "";
	if (f_p) {
		save_prec_file = args::get(f_p);
	}

	// the functions that we are using
	vector<double> xval
	= {0.07768174089481361, 0.4208838472612919,  0.9473139111796449,  0.5089692183323905,
	   0.9570464247595821,  0.15905126169737582, 0.11849894156700524, 0.40999270280625555,
	   0.255893751363913,   0.4596839214805539,  0.48802004339584437, 0.13220719550284032,
	   0.15010376352950006, 0.6220357288995026,  0.5745530849579297,  0.19647252248155644,
	   0.22117387776294362, 0.20684448971820446, 0.24214199341522935, 0.7996623718773789};

	vector<double> yval
	= {0.2626206022801887,  0.5901019591033256, 0.22507728085248468, 0.6589237056448219,
	   0.16792176055765118, 0.7467770932376991, 0.586093747878745,   0.20084730007838192,
	   0.2305544870010452,  0.6032035084811954, 0.7616342948448442,  0.07395053135567398,
	   0.3145108689728848,  0.7994789501256063, 0.31702627032281594, 0.9410740272588111,
	   0.4335486185679469,  0.8896958913102787, 0.8237038487836723,  0.5312785837478274};

	vector<double> alpha = {
	57.12450253383478, 99.4615147296493,  55.26433103539979, 92.3500405527951,  87.67509863441407,
	76.12508998928006, 73.05691678209054, 84.58106997693446, 77.51740813583851, 89.11178685086824,
	79.61010770157222, 77.8956019831367,  76.6795489632914,  71.36011508633008, 99.5339946457969,
	59.75741874591391, 77.37272321552643, 91.75611892120241, 85.7168895671343,  68.55489734818805};

	function<double(double,double)> ffun;
    function<double(double,double)> gfun;
	function<double(double,double)> nfunx;
    function<double(double,double)> nfuny;

	if (f_gauss) {
		ffun = [xval, yval, alpha](double x, double y) {
			double retval = 0;
			for (int i = 0; i < 20; i++) {
				double xv = xval[i];
				double yv = yval[i];
				double a  = alpha[i];
				double r2 = (xv - x) * (xv - x) + (yv - y) * (yv - y);
				retval += 4 * a * (a * r2 - 1) * exp(-a * r2);
			}
			return retval;
		};

		gfun = [xval, yval, alpha](double x, double y) {
			double retval = 0;
			for (int i = 0; i < 20; i++) {
				double xv = xval[i];
				double yv = yval[i];
				double a  = alpha[i];
				double r2 = (xv - x) * (xv - x) + (yv - y) * (yv - y);
				retval += exp(-a * r2);
			}
			return retval;
		};

		nfunx = [xval, yval, alpha](double x, double y) {
			double retval = 0;
			for (int i = 0; i < 20; i++) {
				double xv = xval[i];
				double yv = yval[i];
				double a  = alpha[i];
				double r2 = (xv - x) * (xv - x) + (yv - y) * (yv - y);
				retval += 2 * a * (xv - x) * exp(-a * r2);
			}
			return retval;
		};

		nfuny = [xval, yval, alpha](double x, double y) {
			double retval = 0;
			for (int i = 0; i < 20; i++) {
				double xv = xval[i];
				double yv = yval[i];
				double a  = alpha[i];
				double r2 = (xv - x) * (xv - x) + (yv - y) * (yv - y);
				retval += 2 * a * (yv - y) * exp(-a * r2);
			}
			return retval;
		};
	} else {
		ffun
		= [](double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); };
		gfun = [](double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); };
		nfunx = [](double x, double y) { return M_PI * cos(M_PI * x) * cos(2 * M_PI * y); };
		nfuny = [](double x, double y) { return -2 * M_PI * sin(M_PI * x) * sin(2 * M_PI * y); };
	}

	MPI_Barrier(MPI_COMM_WORLD);
	steady_clock::time_point domain_start = steady_clock::now();

	int              total_domains = num_domains_x * num_domains_y;
	DomainCollection dc(total_domains * my_global_rank / num_procs,
	                    total_domains * (my_global_rank + 1) / num_procs - 1, nx, ny, num_domains_x,
	                    num_domains_y, h_x, h_y, comm);

	if (f_neumann) {
		dc.initNeumann(ffun, gfun, nfunx, nfuny);
	} else {
		dc.initDirichlet(ffun, gfun);
	}

	MPI_Barrier(MPI_COMM_WORLD);
    steady_clock::time_point domain_stop = steady_clock::now();
	duration<double> domain_time = domain_stop - domain_start;

	if (my_global_rank == 0) cout << "Domain Initialization Time: " << domain_time.count() << "\n";

	// Create a map that will be used in the iterative solver
	RCP<map_type> matrix_map = dc.matrix_map;

	// Create the gamma and diff vectors
	RCP<vector_type> gamma = rcp(new vector_type(matrix_map, 1));
	RCP<vector_type> diff  = rcp(new vector_type(matrix_map, 1));

	if (num_domains_x * num_domains_y != 1) {
		// do iterative solve

		// Get the b vector
		RCP<vector_type> b = rcp(new vector_type(matrix_map, 1));
		dc.solveWithInterface(*gamma, *b);

		if (save_rhs_file != "") {
			Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_rhs_file, b, "", "");
		}

		RCP<Tpetra::Operator<>> op;

		if (f_wrapper) {
			// Create a function wrapper
			op = rcp(new FuncWrap(b, &dc));
		} else {
			// Form the matrix
			MPI_Barrier(MPI_COMM_WORLD);
			steady_clock::time_point form_start = steady_clock::now();

			RCP<matrix_type> A = dc.formMatrix(matrix_map);

			MPI_Barrier(MPI_COMM_WORLD);
			duration<double> form_time = steady_clock::now() - form_start;

			if (my_global_rank == 0) cout << "Matrix Formation Time: " << form_time.count() << "\n";

			if (save_matrix_file != "") {
				MPI_Barrier(MPI_COMM_WORLD);
				steady_clock::time_point write_start = steady_clock::now();

				Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(save_matrix_file, A, "",
				                                                           "");

				MPI_Barrier(MPI_COMM_WORLD);
				duration<double> write_time = steady_clock::now() - write_start;
				if (my_global_rank == 0)
					cout << "Time to write matix to file: " << write_time.count() << "\n";
			}

			op = A;
		}

		// Create linear problem for the Belos solver
		Belos::LinearProblem<double, vector_type, Tpetra::Operator<>> problem(op, gamma, b);

		/*toif (f_neumann) {
			RCP<Tpetra::Operator<>> P;
			P = rcp(new ZeroAvg(matrix_map));
			problem.setLeftPrec(P);

		} else*/ if (f_prec) {
			// form preconditioner
			MPI_Barrier(MPI_COMM_WORLD);
			steady_clock::time_point prec_start = steady_clock::now();

			RCP<matrix_type> P = dc.formInvDiag(matrix_map);
            problem.setLeftPrec(P);

			MPI_Barrier(MPI_COMM_WORLD);
			duration<double> prec_time = steady_clock::now() - prec_start;

			if (my_global_rank == 0)
				cout << "Preconditioner Formation Time: " << prec_time.count() << "\n";

			if (save_prec_file != "") {
				MPI_Barrier(MPI_COMM_WORLD);
				steady_clock::time_point write_start = steady_clock::now();

				Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(save_prec_file, P, "",
				                                                           "");

				MPI_Barrier(MPI_COMM_WORLD);
				duration<double> write_time = steady_clock::now() - write_start;
				if (my_global_rank == 0)
					cout << "Time to write preconditioner to file: " << write_time.count() << "\n";
			}
		}

		problem.setProblem();

		MPI_Barrier(MPI_COMM_WORLD);
		steady_clock::time_point iter_start = steady_clock::now();
		// Set the parameters
		Teuchos::ParameterList belosList;
		belosList.set("Block Size", 1);
		belosList.set("Maximum Iterations", 1000);
		belosList.set("Convergence Tolerance", 1e-10);
		int verbosity = Belos::Errors + Belos::StatusTestDetails + Belos::Warnings
		                + Belos::TimingDetails + Belos::Debug;
		belosList.set("Verbosity", verbosity);

		// Create solver and solve
		RCP<Belos::SolverManager<double, vector_type, Tpetra::Operator<>>> solver
		= rcp(new Belos::BlockCGSolMgr<double, vector_type, Tpetra::Operator<>>(
		rcp(&problem, false), rcp(&belosList, false)));
		solver->solve();
		MPI_Barrier(MPI_COMM_WORLD);
		duration<double> iter_time = steady_clock::now() - iter_start;
		if (my_global_rank == 0) std::cout << "CG Time: " << iter_time.count() << "\n";
	}

	// Do one last solve
	MPI_Barrier(MPI_COMM_WORLD);
	steady_clock::time_point solve_start = steady_clock::now();

	dc.solveWithInterface(*gamma, *diff);

	MPI_Barrier(MPI_COMM_WORLD);
	duration<double> solve_time = steady_clock::now() - solve_start;

	if (my_global_rank == 0)
		std::cout << "Time to solve with given set of gammas: " << solve_time.count() << "\n";

	// Calcuate error
	RCP<map_type>    err_map = rcp(new map_type(-1, 1, 0, comm));
	Tpetra::Vector<> exact_norm(err_map);
	Tpetra::Vector<> diff_norm(err_map);

	if (f_neumann) {
		double usum = dc.uSum();
		double uavg;
		Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &usum, &uavg);
		uavg /= total_cells;
		double esum = dc.exactSum();
		double eavg;
		Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &esum, &eavg);
		eavg /= total_cells;

		if (my_global_rank == 0) {
			cerr << "Average of computed solution: " << uavg << "\n";
			cerr << "Average of exact solution: " << eavg << "\n";
		}

		exact_norm.getDataNonConst()[0] = dc.exactNorm(eavg);
		diff_norm.getDataNonConst()[0]  = dc.diffNorm(uavg, eavg);
	} else {
		exact_norm.getDataNonConst()[0] = dc.exactNorm();
		diff_norm.getDataNonConst()[0]  = dc.diffNorm();
	}
	double global_diff_norm;
	double global_exact_norm;
	global_diff_norm  = diff_norm.norm2();
	global_exact_norm = exact_norm.norm2();

	MPI_Barrier(MPI_COMM_WORLD);
	steady_clock::time_point total_stop = steady_clock::now();
	duration<double>         total_time = total_stop - domain_start;

	if (my_global_rank == 0) {
		std::cout << "Total run time: " << total_time.count() << "\n";
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error: " << global_diff_norm / global_exact_norm << "\n";
	}

	MPI_Finalize();
	return 0;
}
