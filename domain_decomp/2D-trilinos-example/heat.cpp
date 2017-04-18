#include "BelosCGIter.hpp"
#include "BelosOutputManager.hpp"
#include "DomainSignatureCollection.h"
#include "FunctionWrapper.h"
#include "MyTypeDefs.h"
#include "ZeroSum.h"
#include "args.h"
//#include <Amesos2.hpp>
//#include <Amesos2_Version.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
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
//#include <mpi.h>
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

	Teuchos::GlobalMPISession global(&argc, &argv, nullptr);
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

	int num_procs = comm->getSize();

	int my_global_rank = comm->getRank();

	// parse input
	args::ArgumentParser  parser("");
	args::HelpFlag        help(parser, "help", "Display this help menu", {'h', "help"});

	args::ValueFlag<int> f_n(parser, "n", "number of cells in the x direction, in each domain",
	                          {'n'});
	args::ValueFlag<string> f_mesh(parser, "file_name", "read in a mesh", {"mesh"});
	args::ValueFlag<int> f_square(
	parser, "num_domains", "create a num_domains x num_domains square of grids", {"square"});
	args::ValueFlag<int> f_amr(parser, "num_domains", "create a num_domains x num_domains square "
	                                                   "of grids, and a num_domains*2 x "
	                                                   "num_domains*2 refined square next to it",
	                            {"amr"});
	args::Flag           f_outclaw(parser, "outclaw", "output amrclaw ascii file", {"outclaw"});
	args::ValueFlag<int> f_l(parser, "n", "run the program n times and print out the average",
	                         {'l'});
	args::ValueFlag<string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                            {'m'});
	args::ValueFlag<string> f_s(parser, "solution filename", "the file to write the solution to",
	                            {'s'});
	args::ValueFlag<string> f_resid(parser, "residual filename",
	                                "the file to write the residual to", {"residual"});
	args::ValueFlag<string> f_error(parser, "error filename",
	                                "the file to write the error to", {"error"});
	args::ValueFlag<string> f_r(parser, "rhs filename", "the file to write the rhs vector to",
	                            {'r'});
	args::ValueFlag<string> f_g(parser, "gamma filename", "the file to write the gamma vector to",
	                            {'g'});
	args::ValueFlag<string> f_p(parser, "preconditioner filename",
	                            "the file to write the preconditioner to", {'p'});
	args::ValueFlag<double> f_t(
	parser, "tolerance", "set the tolerance of the iterative solver (default is 1e-10)", {'t'});
	args::Flag f_wrapper(parser, "wrapper", "use a function wrapper", {"wrap"});
	args::Flag f_gauss(parser, "gauss", "solve gaussian function", {"gauss"});
	args::Flag f_prec(parser, "prec", "use block diagonal preconditioner", {"prec"});
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {"neumann"});
	args::Flag f_cg(parser, "gmres", "use CG for iterative solver", {"cg"});
	args::Flag f_gmres(parser, "gmres", "use GMRES for iterative solver", {"gmres"});
	args::Flag f_rgmres(parser, "rgmres", "use GCRO-DR (Recycling GMRES) for iterative solver",
	                    {"rgmres"});
	args::Flag f_bicg(parser, "gmres", "use BiCGStab for iterative solver", {"bicg"});
	args::Flag f_nozero(parser, "nozero", "don't make the average of vector zero in CG solver",
	                    {"nozero"});
	args::Flag f_lu(parser, "lu", "use LU decomposition", {"lu"});
	args::Flag f_ilu(parser, "ilu", "use incomplete LU preconditioner", {"ilu"});
	args::Flag f_iter(parser, "iterative", "use iterative method", {"iterative"});

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
	DomainSignatureCollection dsc;
	if (f_mesh) {
		string d = args::get(f_mesh);
		dsc      = DomainSignatureCollection(d, comm->getRank());
	} else if (f_amr) {
		int d = args::get(f_amr);
		dsc   = DomainSignatureCollection(d, d, comm->getRank(), true);
	} else {
        int d = args::get(f_square);
		dsc = DomainSignatureCollection(d, d, comm->getRank());
	}
	// Set the number of discretization points in the x and y direction.
	int    nx            = args::get(f_n);
	int    ny            = args::get(f_n);
	int    total_cells   = dsc.num_global_domains*nx*ny;
	if (f_amr) {
		total_cells *= 5;
	}

	if (dsc.num_global_domains < num_procs) {
		std::cerr << "number of domains must be greater than or equal to the number of processes\n";
		return 1;
	}
	// partition domains
	if (num_procs > 1) {
		dsc.zoltanBalance();
	}

	double tol = 1e-10;
	if (f_t) {
		tol = args::get(f_t);
	}

	int loop_count = 1;
	if (f_l) {
		loop_count = args::get(f_l);
	}

	string save_matrix_file = "";
	if (f_m) {
		save_matrix_file = args::get(f_m);
	}

	string save_solution_file = "";
	if (f_s) {
		save_solution_file = args::get(f_s);
	}

	string save_residual_file = "";
	if (f_resid) {
		save_residual_file = args::get(f_resid);
	}
	string save_error_file = "";
	if (f_error) {
		save_error_file = args::get(f_error);
	}
	string save_rhs_file = "";
	if (f_r) {
		save_rhs_file = args::get(f_r);
	}
	string save_gamma_file = "";
	if (f_g) {
		save_gamma_file = args::get(f_g);
	}

	string save_prec_file = "";
	if (f_p) {
		save_prec_file = args::get(f_p);
	}

    
	// the functions that we are using
	function<double(double, double)> ffun;
	function<double(double, double)> gfun;
	function<double(double, double)> nfunx;
	function<double(double, double)> nfuny;

	if (f_gauss) {
		gfun = [](double x, double y) { return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y)); };
		ffun = [](double x, double y) {
			return 100 * M_PI * M_PI * (pow(sin(10 * M_PI * x), 2) - cos(10 * M_PI * x))
			       * exp(cos(10 * M_PI * x))
			       + 121 * M_PI * M_PI * (cos(11 * M_PI * y) - pow(sin(11 * M_PI * y), 2))
			         * exp(cos(11 * M_PI * y));
		};
		nfunx = [](double x, double y) {
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};

		nfuny = [](double x, double y) {
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
	} else {
		ffun = [](double x, double y) {
			x -= .3;
			return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
		};
		gfun = [](double x, double y) {
			x -= .3;
			return sin(M_PI * y) * cos(2 * M_PI * x);
		};
		nfunx = [](double x, double y) {
			x -= .3;
			return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
		};
		nfuny = [](double x, double y) {
			x -= .3;
			return M_PI * cos(M_PI * y) * cos(2 * M_PI * x);
		};
	}

	valarray<double> times(loop_count);
	for (int loop = 0; loop < loop_count; loop++) {
		comm->barrier();
		
		steady_clock::time_point domain_start = steady_clock::now();

		DomainCollection dc(dsc, nx, comm);

		ZeroSum zs;
		if (f_neumann) {
			if (!f_nozero) {
				zs.setTrue();
			}
			dc.initNeumann(ffun, gfun, nfunx, nfuny, f_amr);
		} else {
			dc.initDirichlet(ffun, gfun);
            dc.amr = f_amr;
		}

		comm->barrier();
		steady_clock::time_point domain_stop = steady_clock::now();
		duration<double>         domain_time = domain_stop - domain_start;

		if (my_global_rank == 0)
			cout << "Domain Initialization Time: " << domain_time.count() << "\n";

		// Create a map that will be used in the iterative solver
		RCP<map_type> matrix_map = dc.matrix_map;

		// Create the gamma and diff vectors
		RCP<vector_type> gamma = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type> r     = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type> x     = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type> d     = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type> diff  = rcp(new vector_type(matrix_map, 1));
		RCP<RBMatrix>    RBA;
		RCP<Tpetra::Operator<>> op;

		// Create linear problem for the Belos solver
		RCP<Belos::LinearProblem<double, vector_type, Tpetra::Operator<>>> problem;
		RCP<Belos::SolverManager<double, vector_type, Tpetra::Operator<>>> solver;
		if (dsc.num_global_domains != 1) {
			// do iterative solve

			// Get the b vector
			RCP<vector_type> b = rcp(new vector_type(matrix_map, 1));
			dc.solveWithInterface(*gamma, *b);

			if (save_rhs_file != "") {
				Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_rhs_file, b, "", "");
			}

			if (f_wrapper) {
				// Create a function wrapper
				op = rcp(new FuncWrap(b, &dc));
			} else { // if (f_rbmatrix || f_lu) {
				// Form the matrix
				comm->barrier();
				steady_clock::time_point form_start = steady_clock::now();

				RBA = dc.formRBMatrix(matrix_map);

				comm->barrier();
				duration<double> form_time = steady_clock::now() - form_start;

				if (my_global_rank == 0)
					cout << "Matrix Formation Time: " << form_time.count() << "\n";

				if (save_matrix_file != "") {
					comm->barrier();
					steady_clock::time_point write_start = steady_clock::now();

					ofstream out_file(save_matrix_file);
					out_file << *RBA;
					out_file.close();

					comm->barrier();
					duration<double> write_time = steady_clock::now() - write_start;
					if (my_global_rank == 0)
						cout << "Time to write matrix to file: " << write_time.count() << "\n";
				}
				op = RBA;
			}

			problem
			= rcp(new Belos::LinearProblem<double, vector_type, Tpetra::Operator<>>(op, gamma, b));

			if (f_prec || f_ilu) {
				if (f_ilu) {
					comm->barrier();
					steady_clock::time_point prec_start = steady_clock::now();

					RCP<RBMatrix> L, U;
					RBMatrix      Copy = *RBA;
					Copy.ilu(L, U);
					RCP<LUSolver> solver = rcp(new LUSolver(L, U));
					problem->setLeftPrec(solver);

					comm->barrier();
					duration<double> prec_time = steady_clock::now() - prec_start;

					if (my_global_rank == 0)
						cout << "Preconditioner Formation Time: " << prec_time.count() << "\n";

					// ofstream out_file("L.mm");
					// out_file << *L;
					// out_file.close();
					// out_file = ofstream("U.mm");
					// out_file << *U;
					// out_file.close();
				} else {
					comm->barrier();
					steady_clock::time_point prec_start = steady_clock::now();

					RCP<RBMatrix> P = RBA->invBlockDiag();
					problem->setRightPrec(P);

					comm->barrier();
					duration<double> prec_time = steady_clock::now() - prec_start;

					if (my_global_rank == 0)
						cout << "Preconditioner Formation Time: " << prec_time.count() << "\n";

					if (save_prec_file != "") {
						comm->barrier();
						steady_clock::time_point write_start = steady_clock::now();

						ofstream out_file(save_prec_file);
						out_file << *P;
						out_file.close();

						comm->barrier();
						duration<double> write_time = steady_clock::now() - write_start;
						if (my_global_rank == 0)
							cout << "Time to write preconditioner to file: " << write_time.count()
							     << "\n";
					}
				}
			}

			if (f_lu) {
				comm->barrier();
				steady_clock::time_point lu_start = steady_clock::now();

				RCP<RBMatrix> L, U;
				RBA->lu(L,U);
				LUSolver solver(L,U);
				solver.apply(*b, *gamma);

				comm->barrier();
				duration<double> lu_time = steady_clock::now() - lu_start;
				if (my_global_rank == 0) std::cout << "LU Time: " << lu_time.count() << "\n";

				// ofstream out_file("L.mm");
				// out_file << *L;
				// out_file.close();
				// out_file = ofstream("U.mm");
				// out_file << *U;
				// out_file.close();

			} else {
				problem->setProblem();

				comm->barrier();
				steady_clock::time_point iter_start = steady_clock::now();
				// Set the parameters
				Teuchos::ParameterList belosList;
				belosList.set("Block Size", 1);
				belosList.set("Maximum Iterations", 5000);
				belosList.set("Convergence Tolerance", tol);
				int verbosity = Belos::Errors + Belos::StatusTestDetails + Belos::Warnings
				                + Belos::TimingDetails + Belos::Debug;
				belosList.set("Verbosity", verbosity);

				// Create solver and solve
                if(f_rgmres){
					solver
					= rcp(new Belos::GCRODRSolMgr<double, vector_type, Tpetra::Operator<>>(
					problem, rcp(&belosList, false)));
				} else if (f_cg) {
					solver = rcp(new Belos::BlockCGSolMgr<double, vector_type, Tpetra::Operator<>>(
					problem, rcp(&belosList, false)));
				} else if (f_bicg) {
					solver = rcp(new Belos::BiCGStabSolMgr<double, vector_type, Tpetra::Operator<>>(
					problem, rcp(&belosList, false)));
				} else {
					solver
					= rcp(new Belos::BlockGmresSolMgr<double, vector_type, Tpetra::Operator<>>(
					problem, rcp(&belosList, false)));
				}
				solver->solve();

				comm->barrier();
				duration<double> iter_time = steady_clock::now() - iter_start;
				if (my_global_rank == 0)
					std::cout << "Gamma Solve Time: " << iter_time.count() << "\n";
			}
		}

		// Do one last solve
		comm->barrier();
		steady_clock::time_point solve_start = steady_clock::now();

		dc.solveWithInterface(*gamma, *diff);

		comm->barrier();
		duration<double> solve_time = steady_clock::now() - solve_start;

		if (my_global_rank == 0)
			std::cout << "Time to solve with given set of gammas: " << solve_time.count() << "\n";

        if(f_iter){
            dc.residual();
            dc.swapResidSol();
			dc.solveWithInterface(*x, *r);
            solver->reset(Belos::ResetType::Problem);
            if(f_wrapper){
				((FuncWrap*)op.getRawPtr())->setB(r);
            }
			problem->setProblem(x, r);
			solver->setProblem(problem);
            //cerr<<"solve with 0!!!!!!!!!!!!!!!!!!!!" <<endl;
			solver->solve();
			dc.solveWithInterface(*x, *d);
            dc.sumResidIntoSol();
		}

		// Calcuate error
		RCP<map_type>    err_map = rcp(new map_type(-1, 1, 0, comm));
		Tpetra::Vector<> exact_norm(err_map);
		Tpetra::Vector<> diff_norm(err_map);

		if (f_neumann) {
			double usum = dc.uSum();
			double uavg;
			Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &usum, &uavg);
			uavg /= dc.getGlobalNumCells();
			double esum = dc.exactSum();
			double eavg;
			Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &esum, &eavg);
			eavg /= dc.getGlobalNumCells();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << "\n";
				cout << "Average of exact solution: " << eavg << "\n";
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

		comm->barrier();
		steady_clock::time_point total_stop = steady_clock::now();
		duration<double>         total_time = total_stop - domain_start;

		double residual = dc.residual();
		double fnorm    = dc.fNorm();
		if (my_global_rank == 0) {
			std::cout << "Total run time: " << total_time.count() << "\n";
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << global_diff_norm / global_exact_norm << "\n";
			std::cout << "Residual: " << residual / fnorm << "\n";
			cout << std::fixed;
			cout.precision(2);
		}
		times[loop] = total_time.count();
		if (save_solution_file != "") {
			ofstream out_file(save_solution_file);
			dc.outputSolution(out_file);
			out_file.close();
			if (f_amr) {
				ofstream out_file(save_solution_file + ".amr");
				dc.outputSolutionRefined(out_file);
				out_file.close();
			}
		}
		if (save_residual_file != "") {
			ofstream out_file(save_residual_file);
			dc.outputResidual(out_file);
			out_file.close();
			if (f_amr) {
				ofstream out_file(save_residual_file + ".amr");
				dc.outputResidualRefined(out_file);
				out_file.close();
			}
		}
		if (save_error_file != "") {
			ofstream out_file(save_error_file);
			dc.outputError(out_file);
			out_file.close();
			if (f_amr) {
				ofstream out_file(save_error_file + ".amr");
				dc.outputErrorRefined(out_file);
				out_file.close();
			}
		}
		if (save_gamma_file != "") {
			Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_gamma_file, gamma, "",
			                                                          "");
		}
		if (f_outclaw) {
			dc.outputClaw();
		}
	}

	if (loop_count > 1 && my_global_rank == 0) {
		cout << std::fixed;
		cout.precision(2);
		std::cout << "Times: ";
		for (double t : times) {
			cout << t << " ";
		}
		cout << "\n";
		cout << "Average: " << times.sum() / times.size() << "\n";
	}

	return 0;
}
