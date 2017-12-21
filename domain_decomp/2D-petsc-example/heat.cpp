#include "DomainCollection.h"
#include "FivePtPatchOperator.h"
#include "FourthInterpolator.h"
//#include "FunctionWrapper.h"
#include "Init.h"
#include "MatrixHelper.h"
#include "PW.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "PatchSolvers/FishpackPatchSolver.h"
#include "QuadInterpolator.h"
#include "SchurHelper.h"
#include "Writers/ClawWriter.h"
#include "Writers/MMWriter.h"
#ifdef ENABLE_AMGX
#include "AmgxWrapper.h"
#endif
#ifdef HAVE_VTK
#include "Writers/VtkWriter.h"
#endif
#include "Timer.h"
#include "args.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <sstream>
#include <string>
#include <string>

// =========== //
// main driver //
// =========== //

using namespace std;

int main(int argc, char *argv[])
{
	int     petsc_argc;
	int *   petsc_argc_ptr = nullptr;
	char ** petsc_argv;
	char ***petsc_argv_ptr = nullptr;
	string  delim          = "::";
	for (int i = 0; i < argc; i++) {
		if (argv[i] == delim) {
			int tmp        = argc;
			argc           = i;
			petsc_argc     = tmp - i;
			petsc_argc_ptr = &petsc_argc;
			petsc_argv     = &argv[i];
			petsc_argv_ptr = &petsc_argv;
		}
	}

	// use :: to delimit petsc options at end
	PetscInitialize(petsc_argc_ptr, petsc_argv_ptr, nullptr, nullptr);

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

	// patch solvers
	args::Flag f_fish(parser, "fishpack", "use fishpack as the patch solver", {"fishpack"});
#ifdef __NVCC__
	args::Flag f_cufft(parser, "cufft", "use CuFFT as the patch solver", {"cufft"});
#endif

	// third-party preconditioners
	args::Flag f_amgx(parser, "amgx", "solve schur compliment system with amgx", {"amgx"});

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

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
		dc       = DomainCollection(d);
	} else if (f_amr) {
		int d = args::get(f_amr);
		dc    = DomainCollection(d, d, true);
	} else {
		int d = args::get(f_square);
		dc    = DomainCollection(d, d);
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
	double tol = 1e-12;
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
	shared_ptr<PatchSolver> p_solver(new FftwPatchSolver(dc, lambda));

	// patch operator
	shared_ptr<PatchOperator> p_operator(new FivePtPatchOperator());

	// interface interpolator
	shared_ptr<Interpolator> p_interp(new QuadInterpolator());

	Tools::Timer timer;
	timer.start("Domain Initialization");

	SchurHelper  sch(dc, p_solver, p_operator, p_interp);
	MatrixHelper mh(dc);

	PW<Vec> u_next = dc.getNewDomainVec();
	PW<Vec> u      = dc.getNewDomainVec();
	PW<Vec> u_prev = dc.getNewDomainVec();
	PW<Vec> u_tmp;
	PW<Vec> f = dc.getNewDomainVec();

	double t = 0;
	if (f_bdf2) {
		Init::fillSolution(dc, u_prev, efun, 0);
		Init::fillSolution(dc, u, efun, dt);
		t = dt;
	} else {
		Init::fillSolution(dc, u, efun, 0);
	}

	timer.stop("Domain Initialization");

	// Create the gamma and diff vectors
	PW<Vec> gamma   = dc.getNewSchurVec();
	PW<Vec> gamma_e = dc.getNewSchurVec();
	PW<Vec> zeros   = dc.getNewSchurVec();
	PW<Vec> diff    = dc.getNewSchurVec();
	PW<Vec> b       = dc.getNewSchurVec();
	PW<Mat> A;

	// Create linear solver
	PW<KSP> solver;
	KSPCreate(MPI_COMM_WORLD, &solver);
	KSPSetFromOptions(solver);

	PW<Vec> exact = dc.getNewDomainVec();
	PW<Vec> error = dc.getNewDomainVec();
	PW<Vec> resid = dc.getNewDomainVec();
	if ((f_schur && dc.num_global_domains != 1) || !f_schur) {
		// do iterative solve

		///////////////////
		// setup start
		///////////////////
		timer.start("Linear System Setup");

		timer.start("Matrix Formation");

		if (f_schur) {
			A = sch.formCRSMatrix();
			if (f_wrapper) {
				// Create a function wrapper
				// op = rcp(new FuncWrap(&b, &u_next, &f, &sch));
			}
		} else {
			A = mh.formCRSMatrix(lambda);
		}

		timer.stop("Matrix Formation");
		if (f_m) {
			PetscViewer viewer;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, args::get(f_m).c_str(), FILE_MODE_WRITE,
			                      &viewer);
			MatView(A, viewer);
			PetscViewerDestroy(&viewer);
		}

		KSPSetTolerances(solver, tol, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
		KSPSetOperators(solver, A, A);
		KSPSetUp(solver);
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

	int i = 0;

	for (; t < tend; t += dt) {
#ifdef HAVE_VTK
		if (f_outvtk) {
			if (t >= i * out_every) {
				std::ostringstream ss;
				ss << args::get(f_outvtk) << setfill('0') << setw(5) << i;
				VtkWriter writer(dc, ss.str());
				writer.add(u, "Solution");
				Init::fillSolution(dc, exact, efun, t - dt);
				writer.add(exact, "Exact");
				VecAXPBYPCZ(error, -1.0, 1.0, 0.0, exact, u);
				writer.add(error, "Error");
				if (f_schur) {
					sch.applyWithInterface(u, gamma, resid);
				} else {
					MatMult(A, u, resid);
				}
				VecAYPX(resid, -1.0, f);
				writer.add(resid, "Residual");
				writer.write();
				i++;
			}
		}
#endif
		if (f_bdf2) {
			VecAXPBYPCZ(f, -1.0 / 3.0, 4.0 / 3.0, 0.0, u_prev, u);
		} else {
			VecScale(f, 0);
			VecAXPY(f, 1.0, u);
		}
		VecScale(f, lambda);

		if ((f_schur && dc.num_global_domains != 1) || !f_schur) {
			timer.start("Linear Solve");
			// Get the b vector
			if (f_schur) {
				sch.solveWithInterface(f, u_next, zeros, b);
				VecScale(b, -1.0);
				if (f_r) {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, args::get(f_r).c_str(), FILE_MODE_WRITE,
					                      &viewer);
					VecView(b, viewer);
					PetscViewerDestroy(&viewer);
				}
			}
			// solve
			if (f_schur) {
				KSPSolve(solver, b, gamma);
			} else {
				KSPSolve(solver, f, u_next);
			}
			timer.stop("Linear Solve");
		}

		// Do one last solve
		if (f_schur) {
			timer.start("Patch Solve");

			// gamma->scale(lambda);
			sch.solveWithInterface(f, u_next, gamma, diff);
			timer.stop("Patch Solve");
		}
		u_tmp  = u_prev;
		u_prev = u;
		u      = u_next;
		u_next = u_tmp;
		// use previous solution as initial guess
		VecAXPY(u_next, 1.0, u);
	}
	///////////////////
	// solve end
	///////////////////
	timer.stop("Complete Solve");
	// error
	Init::fillSolution(dc, exact, efun, t);
	VecAXPBYPCZ(error, -1.0, 1.0, 0.0, exact, u);

	double error_norm, error_norm_inf, exact_norm;
	VecNorm(error, NORM_2, &error_norm);
	VecNorm(error, NORM_INFINITY, &error_norm_inf);
	VecNorm(exact, NORM_2, &exact_norm);

	if (f_schur) {
		sch.applyWithInterface(u, gamma, resid);
	} else {
		MatMult(A, u, resid);
	}
	VecAYPX(resid, -1.0, f);

	double resid_norm, resid_norm_inf, f_norm;
	VecNorm(resid, NORM_2, &resid_norm);
	VecNorm(resid, NORM_INFINITY, &resid_norm_inf);
	VecNorm(f, NORM_2, &f_norm);

	if (my_global_rank == 0) {
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error    (2 norm):   " << error_norm / exact_norm << endl;
		std::cout << "Error    (inf norm): " << error_norm_inf << endl;
		std::cout << "Residual (2 norm):   " << resid_norm / f_norm << endl;
		std::cout << "Residual (inf norm): " << resid_norm_inf << endl;
		cout.unsetf(std::ios_base::floatfield);
	}
	if (f_g) {
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, args::get(f_g).c_str(), &viewer);
		VecView(gamma, viewer);
	}
	if (f_outclaw) {
		ClawWriter writer(dc);
		writer.write(u, u);
	}
#ifdef HAVE_VTK
	if (f_outvtk) {
		std::ostringstream ss;
		ss << args::get(f_outvtk) << setfill('0') << setw(5) << i;
		VtkWriter writer(dc, ss.str());
		writer.add(u, "Solution");
		writer.add(exact, "Exact");
		writer.add(error, "Error");
		writer.add(resid, "Residual");
		writer.write();
	}
#endif
	cout.unsetf(std::ios_base::floatfield);

	if (my_global_rank == 0) {
		cout << timer;
	}
	return 0;
}
