#include "DomainCollection.h"
#include "Factory.h"
#include "FivePtPatchOperator.h"
#include "FourthInterpolator.h"
#include "FunctionWrapper.h"
#include "Init.h"
#include "MatrixHelper.h"
#include "MyTypeDefs.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "PatchSolvers/FishpackPatchSolver.h"
#include "QuadInterpolator.h"
#include "SchurHelper.h"
#include "Writers/ClawWriter.h"
#include "Writers/MMWriter.h"
#ifdef ENABLE_AMGX
#include "AmgxWrapper.h"
#endif
#ifdef ENABLE_HYPRE
#include "HypreWrapper.h"
#endif
#ifdef HAVE_VTK
#include "Writers/VtkWriter.h"
#endif
#include "Timer.h"
#include "args.h"
#include <memory>
#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosLSQRSolMgr.hpp>
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
#include <Tpetra_Experimental_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
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

	Teuchos::GlobalMPISession     global(&argc, &argv, nullptr);
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

	int num_procs = comm->getSize();

	int my_global_rank = comm->getRank();

	// parse input
	args::ArgumentParser parser("");

	// program options
	args::HelpFlag       help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<int> f_n(parser, "n", "number of cells in the x direction, in each domain",
	                         {'n'});
	args::ValueFlag<int> f_l(parser, "n", "run the program n times and print out the average",
	                         {'l'});
	args::Flag           f_wrapper(parser, "wrapper", "use a function wrapper", {"wrap"});
	args::ValueFlag<int> f_div(parser, "divide", "use iterative method", {"divide"});
	args::Flag f_nozerof(parser, "zerou", "don't  make the rhs match the boundary conditions",
	                     {"nozerof"});
	args::ValueFlag<double> f_t(
	parser, "tolerance", "set the tolerance of the iterative solver (default is 1e-10)", {'t'});
	args::ValueFlag<double> f_dt(parser, "tolerance", "time step size", {"dt"});
	args::ValueFlag<double> f_out_every(parser, "tolerance", "time step size", {"outevery"});
	args::ValueFlag<double> f_tend(parser, "", "stop time", {"tend"});
	args::ValueFlag<double> f_alpha(parser, "", "diffusivity coefficient", {"alpha"});
	args::Flag              f_schur(parser, "", "use schur compliment method for solve", {"schur"});
	args::Flag              f_bdf2(parser, "", "use bdf2 time stepping", {"bdf2"});

	// mesh options
	args::ValueFlag<string> f_mesh(parser, "file_name", "read in a mesh", {"mesh"});
	args::ValueFlag<int>    f_square(parser, "num_domains",
	                              "create a num_domains x num_domains square of grids", {"square"});
	args::ValueFlag<int> f_amr(parser, "num_domains", "create a num_domains x num_domains square "
	                                                  "of grids, and a num_domains*2 x "
	                                                  "num_domains*2 refined square next to it",
	                           {"amr"});

	// output options
	args::Flag f_outclaw(parser, "outclaw", "output amrclaw ascii file", {"outclaw"});
#ifdef HAVE_VTK
	args::ValueFlag<string> f_outvtk(parser, "", "output to vtk format", {"outvtk"});
#endif
	args::ValueFlag<string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                            {'m'});
	args::ValueFlag<string> f_r(parser, "rhs filename", "the file to write the rhs vector to",
	                            {'r'});
	args::ValueFlag<string> f_g(parser, "gamma filename", "the file to write the gamma vector to",
	                            {'g'});

	// preconditioners
	args::Flag f_precj(parser, "prec", "use jacobi preconditioner", {"precj"});
	args::Flag f_prec(parser, "prec", "use block diagonal jacobi preconditioner", {"prec"});
	args::Flag f_precmuelu(parser, "prec", "use MueLu AMG preconditioner", {"muelu"});
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {"neumann"});
	args::Flag f_cg(parser, "gmres", "use CG for iterative solver", {"cg"});
	args::Flag f_gmres(parser, "gmres", "use GMRES for iterative solver", {"gmres"});
	args::Flag f_lsqr(parser, "gmres", "use least squares for iterative solver", {"lsqr"});
	args::Flag f_rgmres(parser, "rgmres", "use GCRO-DR (Recycling GMRES) for iterative solver",
	                    {"rgmres"});
	args::Flag f_bicg(parser, "gmres", "use BiCGStab for iterative solver", {"bicg"});

	// direct solvers
	args::Flag f_iter(parser, "iterative", "use iterative method", {"iterative"});

	// patch solvers
	args::Flag f_fish(parser, "fishpack", "use fishpack as the patch solver", {"fishpack"});
#ifdef __NVCC__
	args::Flag f_cufft(parser, "cufft", "use CuFFT as the patch solver", {"cufft"});
#endif

	// third-party preconditioners
	args::Flag f_amgx(parser, "amgx", "solve schur compliment system with amgx", {"amgx"});
	args::Flag f_hypre(parser, "hypre", "solve schur compliment system with hypre", {"hypre"});

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

	DomainCollection dc;
	if (f_mesh) {
		string d = args::get(f_mesh);
		dc       = DomainCollection(comm, d, comm->getRank());
	} else if (f_amr) {
		int d = args::get(f_amr);
		dc    = DomainCollection(comm, d, d, comm->getRank(), true);
	} else {
		int d = args::get(f_square);
		dc    = DomainCollection(comm, d, d, comm->getRank());
	}
	if (f_div) {
		for (int i = 0; i < args::get(f_div); i++) {
			dc.divide();
		}
	}
	// Set the number of discretization points in the x and y direction.
	int nx = args::get(f_n);
	int ny = args::get(f_n);
	dc.n   = nx;
	for (auto &p : dc.domains) {
		p.second.n = nx;
	}
	int total_cells = dc.num_global_domains * nx * ny;
	cerr << "Total cells: " << total_cells << endl;

	if (dc.num_global_domains < num_procs) {
		std::cerr << "number of domains must be greater than or equal to the number of processes\n";
		return 1;
	}
	// partition domains
	if (num_procs > 1) {
		dc.zoltanBalance();
	}

	double dt        = args::get(f_dt);
	double tend      = args::get(f_tend);
	double out_every = 0;
	if (f_out_every) {
		out_every = args::get(f_out_every);
	}
	scalar_type tol = 1e-12;
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

	string save_rhs_file = "";
	if (f_r) {
		save_rhs_file = args::get(f_r);
	}
	string save_gamma_file = "";
	if (f_g) {
		save_gamma_file = args::get(f_g);
	}

	double alpha = 0.5;
	if (f_alpha) {
		alpha = args::get(f_alpha);
	}
	// the functions that we are using

	auto efun = [=](double x, double y, double t) {
		return exp(-2 * M_PI * M_PI * alpha * t) * sin(1 * M_PI * y) * sin(1 * M_PI * x);
	};
	// set the patch solver
	double lambda = -1.0 / (alpha * dt);
	if (f_bdf2) {
		lambda = -3.0 / (2.0 * alpha * dt);
	}
	RCP<PatchSolver> p_solver = rcp(new FftwPatchSolver(dc, lambda));

	// patch operator
	RCP<PatchOperator> p_operator = rcp(new FivePtPatchOperator());

	// interface interpolator
	RCP<Interpolator> p_interp = rcp(new QuadInterpolator());

	Tools::Timer timer;
	timer.start("Domain Initialization");

	SchurHelper  sch(dc, comm, p_solver, p_operator, p_interp);
	MatrixHelper mh(dc, comm);

	IS    domain_is = dc.getDomainIS();
	RCP<map_type>    domain_map = dc.getDomainRowMap();
	RCP<vector_type> u_next     = rcp(new vector_type(domain_map, 1));
	RCP<vector_type> u          = rcp(new vector_type(domain_map, 1));
	RCP<vector_type> u_prev     = rcp(new vector_type(domain_map, 1));
	RCP<vector_type> u_tmp;
	RCP<vector_type> f = rcp(new vector_type(domain_map, 1));

	double t = 0;
	if (f_bdf2) {
		Init::fillSolution(dc, *u_prev, efun, 0);
		Init::fillSolution(dc, *u, efun, dt);
		t = dt;
	} else {
		Init::fillSolution(dc, *u, efun, 0);
	}

	timer.stop("Domain Initialization");

	// Create a map that will be used in the iterative solver
	RCP<map_type>       matrix_map       = dc.getSchurRowMap();
	RCP<const map_type> matrix_map_const = matrix_map;

	// Create the gamma and diff vectors
	RCP<vector_type>                   gamma   = rcp(new vector_type(matrix_map, 1));
	RCP<vector_type>                   gamma_e = rcp(new vector_type(matrix_map, 1));
	RCP<vector_type>                   zeros   = rcp(new vector_type(matrix_map, 1));
	RCP<vector_type>                   diff    = rcp(new vector_type(matrix_map, 1));
	RCP<vector_type>                   b       = rcp(new vector_type(matrix_map, 1));
	RCP<matrix_type>                   A;
	RCP<matrix_type>                   A_full;
	RCP<const Tpetra::RowMatrix<>>     rm;
	RCP<Tpetra::Operator<scalar_type>> op;

	// Create linear problem for the Belos solver
	RCP<Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>> problem;
	RCP<Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>>
	problem_resid;
	RCP<Belos::SolverManager<scalar_type, vector_type, Tpetra::Operator<scalar_type>>> solver;
	Teuchos::ParameterList belosList;

	typedef Ifpack2::Preconditioner<scalar_type> Preconditioner;
	RCP<Preconditioner>                          prec;

#ifdef ENABLE_AMGX
	Teuchos::RCP<AmgxWrapper> amgxsolver;
#endif
#ifdef ENABLE_HYPRE
	Teuchos::RCP<HypreWrapper> hypresolver;
#endif

		RCP<vector_type> exact = rcp(new vector_type(domain_map, 1));
		RCP<vector_type> error = rcp(new vector_type(domain_map, 1));
		RCP<vector_type> resid = rcp(new vector_type(domain_map, 1));
	if ((f_schur && dc.num_global_domains != 1) || !f_schur) {
		// do iterative solve

		///////////////////
		// setup start
		///////////////////
		timer.start("Linear System Setup");

		timer.start("Matrix Formation");

		if (f_schur) {
			A      = sch.formCRSMatrix();
			A_full = mh.formCRSMatrix(lambda);
			rm     = A;
			op     = A;
			if (f_wrapper) {
				// Create a function wrapper
				op = rcp(new FuncWrap(&b, &u_next, &f, &sch));
			}
		} else {
			A      = mh.formCRSMatrix(lambda);
			A_full = A;
			rm     = A;
			op     = A;
		}

		timer.stop("Matrix Formation");

		if (save_matrix_file != "")
			Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(save_matrix_file, A, "", "");

		if (f_schur) {
			problem
			= rcp(new Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>(
			op, gamma, b));

			problem_resid
			= rcp(new Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>(
			op, gamma_e, b));
			if (f_hypre&&f_wrapper) {
					// Create the relaxation.  You could also do this using
					// Ifpack2::Factory (the preconditioner factory) if you like.
					RCP<precond_type> prec = rcp(new Ifpack2::Relaxation<Tpetra::RowMatrix<>>(A));
					// Make the list of relaxation parameters.
					Teuchos::ParameterList params;
					// Do symmetric SOR / Gauss-Seidel.
					params.set("relaxation: type", "Jacobi");
					// Two sweeps (of symmetric SOR / Gauss-Seidel) per apply() call.
					params.set("relaxation: sweeps", 1);
					// ... Set any other parameters you want to set ...

					// Set parameters.
					prec->setParameters(params);
					// Prepare the relaxation instance for use.
					prec->initialize();
					prec->compute();

				problem->setLeftPrec(prec);

			}
		} else {
			problem
			= rcp(new Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>(
			op, u_next, f));
		}
		if (f_amgx) {
#ifdef ENABLE_AMGX
			timer.start("AMGX Setup");
			amgxsolver = rcp(new AmgxWrapper(A, dc, nx));
			timer.stop("AMGX Setup");
#endif
		} else if (f_hypre&&!f_wrapper ) {
#ifdef ENABLE_HYPRE
			timer.start("Hypre Setup");
			hypresolver = rcp(new HypreWrapper(A, dc, nx, tol, true));
			timer.stop("Hypre Setup");
#endif
		} else {
			// Set the parameters
			belosList.set("Block Size", 1);
			belosList.set("Maximum Iterations", 5000);
			belosList.set("Convergence Tolerance", tol);
			belosList.set("Output Frequency", out_every==0?1:(out_every / dt));
			// int verbosity = 0;
			int verbosity = Belos::TimingDetails + Belos::IterationDetails;
			belosList.set("Verbosity", verbosity);
			belosList.set("Rel RHS Err", 0.0);
			belosList.set("Rel Mat Err", 0.0);

			// Create solver and solve
			solver = rcp(
			new Belos::BlockGmresSolMgr<scalar_type, vector_type, Tpetra::Operator<scalar_type>>(
			problem, rcp(&belosList, false)));
		}
		///////////////////
		// setup end
		///////////////////
		timer.stop("Linear System Setup");
	}
	///////////////////
	// solve start
	///////////////////

	// initialize vectors
	timer.start("Complete Solve");

	vector_type ones(domain_map, 1);
	ones.putScalar(1);
	int i = 0;

	for (; t < tend; t += dt) {
#ifdef HAVE_VTK
		if (f_outvtk) {
			if (t >= i * out_every) {
				std::ostringstream ss;
				ss << args::get(f_outvtk) << setfill('0') << setw(5) << i;
				VtkWriter writer(dc, ss.str());
				writer.add(*u, "Solution");
				Init::fillSolution(dc, *exact, efun, t - dt);
				writer.add(*exact, "Exact");
				error->update(-1.0, *exact, 1.0, *u, 0.0);
				writer.add(*error, "Error");
				A_full->apply(*u, *resid);
				resid->update(-1.0, *f, 1.0);
				writer.add(*resid, "Residual");
				writer.write();
				i++;
			}
		}
#endif
		if (f_bdf2) {
			f->update(-1.0 / 3.0, *u_prev, 4.0 / 3.0, *u, 0.0);
		} else {
			f->update(1.0, *u, 0.0);
		}
		f->scale(lambda);

		if ((f_schur && dc.num_global_domains != 1) || !f_schur) {
			timer.start("Linear Solve");
			// Get the b vector
			if (f_schur) {
				sch.solveWithInterface(*f, *u_next, *zeros, *b);
				b->scale(-1.0);
				if (save_rhs_file != "") {
					Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_rhs_file, b, "",
					                                                          "");
				}
			}
			// solve
			if (false) {
#ifdef ENABLE_AMGX
			} else if (f_amgx) {
				if (f_schur) {
					amgxsolver->solve(gamma, b);
				} else {
					amgxsolver->solve(u_next, f);
				}
#endif
#ifdef ENABLE_HYPRE
			} else if (f_hypre&&!f_wrapper) {
				if (f_schur) {
					hypresolver->solve(gamma, b);
				} else {
					hypresolver->solve(u_next, f);
				}
#endif
			} else {
                if(f_schur){
				problem->setProblem(gamma,b);
                }else{
				problem->setProblem(u_next,f);
                }

				solver->solve();
			}
			timer.stop("Linear Solve");
		}

		// Do one last solve
		if (f_schur) {
			timer.start("Patch Solve");

			// gamma->scale(lambda);
			sch.solveWithInterface(*f, *u_next, *gamma, *diff);
			timer.stop("Patch Solve");
			if (f_iter) {
				A_full->apply(*u, *resid);
				resid->update(-1.0, *f, 1.0);
				sch.solveWithInterface(*resid, *error, *gamma_e, *b);
				b->scale(-1);

				problem_resid->setProblem();
				solver = rcp(new Belos::BlockGmresSolMgr<scalar_type, vector_type,
				                                         Tpetra::Operator<scalar_type>>(
				problem_resid, rcp(&belosList, false)));
				solver->solve();
				sch.solveWithInterface(*resid, *error, *gamma_e, *b);
				u_next->update(1, *error, 1);
			}
		}
		u_tmp  = u_prev;
		u_prev = u;
		u      = u_next;
		u_next = u_tmp;
		// use previous solution as initial guess
		u_next->update(1.0, *u, 0.0);
	}
	///////////////////
	// solve end
	///////////////////
	timer.stop("Complete Solve");
	// error
	Init::fillSolution(dc, *exact, efun, t);
	error->update(-1.0, *exact, 1.0, *u, 0.0);
	double error_norm     = error->getVector(0)->norm2();
	double error_norm_inf = error->getVector(0)->normInf();
	double exact_norm     = exact->getVector(0)->norm2();
	A_full->apply(*u, *resid);
	resid->update(-1.0, *f, 1.0);
	double resid_norm     = resid->getVector(0)->norm2();
	double resid_norm_inf = resid->getVector(0)->normInf();
	double f_norm         = f->getVector(0)->norm2();

	if (my_global_rank == 0) {
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error    (2 norm):   " << error_norm / exact_norm << endl;
		std::cout << "Error    (inf norm): " << error_norm_inf << endl;
		std::cout << "Residual (2 norm):   " << resid_norm / f_norm << endl;
		std::cout << "Residual (inf norm): " << resid_norm_inf << endl;
		cout.unsetf(std::ios_base::floatfield);
	}
	if (save_gamma_file != "") {
		Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_gamma_file, gamma, "", "");
	}
	if (f_outclaw) {
		ClawWriter writer(dc);
		writer.write(*u, *u);
	}
#ifdef HAVE_VTK
	if (f_outvtk) {
		std::ostringstream ss;
		ss << args::get(f_outvtk) << setfill('0') << setw(5) << i;
		VtkWriter writer(dc, ss.str());
		writer.add(*u, "Solution");
		writer.add(*exact, "Exact");
		writer.add(*error, "Error");
		writer.add(*resid, "Residual");
		writer.write();
	}
#endif
	cout.unsetf(std::ios_base::floatfield);

	if (my_global_rank == 0) {
		cout << timer;
	}
	return 0;
}
