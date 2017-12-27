#include "DomainCollection.h"
#include "SchurHelper.h"
//#include "FunctionWrapper.h"
#include "FivePtPatchOperator.h"
#include "Init.h"
#include "MatrixHelper.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "QuadInterpolator.h"
#include "Writers/ClawWriter.h"
#include "Writers/MMWriter.h"
#ifdef ENABLE_AMGX
#include "AmgxWrapper.h"
#endif
#ifdef ENABLE_MUELU
#include <MueLuWrapper.h>
#endif
#ifdef ENABLE_MUELU_CUDA
#include <MueLuCudaWrapper.h>
#endif
#ifdef HAVE_VTK
#include "Writers/VtkWriter.h"
#endif
#include "Timer.h"
#include "args.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <string>
#include <unistd.h>

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
	args::Flag f_neumann(parser, "", "use neumann boundary conditions", {"neumann"});
	args::Flag f_noschur(parser, "", "don't use schur compliment method", {"noschur"});

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
	args::ValueFlag<string> f_s(parser, "solution filename", "the file to write the solution to",
	                            {'s'});
	args::ValueFlag<string> f_resid(parser, "residual filename",
	                                "the file to write the residual to", {"residual"});
	args::ValueFlag<string> f_error(parser, "error filename", "the file to write the error to",
	                                {"error"});
	args::ValueFlag<string> f_rhs(parser, "rhs filename", "the file to write the rhs vector to",
	                              {'r'});
	args::ValueFlag<string> f_g(parser, "gamma filename", "the file to write the gamma vector to",
	                            {'g'});
	args::ValueFlag<string> f_read_gamma(parser, "gamma filename",
	                                     "the file to read gamma vector from", {"readgamma"});

	// problem options
	args::Flag f_gauss(parser, "gauss", "solve gaussian function", {"gauss"});
	args::Flag f_zero(parser, "zero", "solve zero function", {"zero"});

// patch solvers
#ifdef __NVCC__
	args::Flag f_cufft(parser, "cufft", "use CuFFT as the patch solver", {"cufft"});
#endif

// third-party preconditioners
#ifdef ENABLE_AMGX
	args::ValueFlag<string> f_amgx(parser, "", "solve schur compliment system with amgx", {"amgx"});
#endif
#ifdef ENABLE_MUELU
	args::ValueFlag<string> f_meulu(parser, "", "solve schur compliment system with muelu",
	                                {"muelu"});
#endif
#ifdef ENABLE_MUELU_CUDA
	args::ValueFlag<string> f_meulucuda(parser, "", "solve schur compliment system with muelu",
	                                    {"muelucuda"});
#endif

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
	// Set the number of discretization points in the x and y direction.
	int n = args::get(f_n);

	double tol = 1e-12;
	if (f_t) {
		tol = args::get(f_t);
	}

	int loop_count = 1;
	if (f_l) {
		loop_count = args::get(f_l);
	}
	/***********
	 * Input parsing done
	 **********/

	///////////////
	// Create Mesh
	///////////////
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
	if (f_neumann) {
		dc.setNeumann();
	}
	dc.n = n;
	for (auto &p : dc.domains) {
		p.second.n = n;
	}

	if (dc.num_global_domains < num_procs) {
		std::cerr << "number of domains must be greater than or equal to the "
		             "number of processes\n";
		return 1;
	}
	// partition domains if running in parallel
	if (num_procs > 1) {
		dc.zoltanBalance();
	}

	// the functions that we are using
	function<double(double, double, double)> ffun;
	function<double(double, double, double)> gfun;
	function<double(double, double, double)> nfunx;
	function<double(double, double, double)> nfuny;
	function<double(double, double, double)> nfunz;

	if (f_zero) {
		ffun  = [](double x, double y, double z) { return 0; };
		gfun  = [](double x, double y, double z) { return 0; };
		nfunx = [](double x, double y, double z) { return 0; };
		nfuny = [](double x, double y, double z) { return 0; };
		nfunz = [](double x, double y, double z) { return 0; };
	} else if (f_gauss) {
		gfun = [](double x, double y, double z) {
			return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y));
		};
		ffun = [](double x, double y, double z) {
			return 100 * M_PI * M_PI * (pow(sin(10 * M_PI * x), 2) - cos(10 * M_PI * x))
			       * exp(cos(10 * M_PI * x))
			       + 121 * M_PI * M_PI * (cos(11 * M_PI * y) - pow(sin(11 * M_PI * y), 2))
			         * exp(cos(11 * M_PI * y));
		};
		nfunx = [](double x, double y, double z) {
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};

		nfuny = [](double x, double y, double z) {
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
		nfunz = [](double x, double y, double z) { return 0; };
	} else {
		ffun = [](double x, double y, double z) {
			return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
		};
		gfun  = [](double x, double y, double z) { return sin(M_PI * y) * cos(2 * M_PI * x); };
		nfunx = [](double x, double y, double z) {
			return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
		};
		nfuny
		= [](double x, double y, double z) { return M_PI * cos(M_PI * y) * cos(2 * M_PI * x); };
		nfunz = [](double x, double y, double z) { return 0; };
	}

	// set the patch solver
	shared_ptr<PatchSolver> p_solver;
	p_solver.reset(new FftwPatchSolver(dc));

	// patch operator
	shared_ptr<PatchOperator> p_operator(new FivePtPatchOperator());

	// interface interpolator
	shared_ptr<Interpolator> p_interp(new QuadInterpolator());

#ifdef ENABLE_AMGX
	AmgxWrapper *amgxsolver = nullptr;
	if (f_amgx) {
		amgxsolver = new AmgxWrapper(args::get(f_amgx));
	}
#endif

#ifdef ENABLE_MUELU_CUDA
	if (f_meulucuda) {
		MueLuCudaWrapper::initialize();
	}
#endif
	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		SchurHelper  sch(dc, p_solver, p_operator, p_interp);
		MatrixHelper mh(dc);

		PW<Vec> u     = dc.getNewDomainVec();
		PW<Vec> exact = dc.getNewDomainVec();
		PW<Vec> f     = dc.getNewDomainVec();

		if (f_neumann) {
			Init::initNeumann(dc, n, f, exact, ffun, gfun, nfunx, nfuny,nfunz);
		} else {
			Init::initDirichlet(dc, n, f, exact, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		// Create the gamma and diff vectors
		PW<Vec> gamma = dc.getNewSchurVec();
		PW<Vec> diff  = dc.getNewSchurVec();
		PW<Vec> b     = dc.getNewSchurVec();
		PW<Mat> A;

		// Create linear problem for the Belos solver
		PW<KSP> solver;
		KSPCreate(MPI_COMM_WORLD, &solver);
		KSPSetFromOptions(solver);

		if (f_neumann && !f_nozerof) {
			double fdiff = dc.integrate(f) / dc.area();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
			VecShift(f, -fdiff);
		}

#ifdef ENABLE_MUELU
		MueLuWrapper *meulusolver = nullptr;
#endif
#ifdef ENABLE_MUELU_CUDA
		MueLuCudaWrapper *meulucudasolver = nullptr;
#endif

		if (f_noschur || dc.num_global_domains != 1) {
			// do iterative solve

			if (!f_noschur) {
				// Get the b vector
				VecScale(gamma, 0);
				sch.solveWithInterface(f, u, gamma, b);
				VecScale(b, -1.0);

				if (f_rhs) {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, args::get(f_rhs).c_str(),
					                      FILE_MODE_WRITE, &viewer);
					VecView(b, viewer);
					PetscViewerDestroy(&viewer);
				}
			}

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			if (f_wrapper) {
				// Create a function wrapper
				// op = rcp(new FuncWrap(b, &dc));
			} else {
				timer.start("Matrix Formation");

				if (f_noschur) {
					A = mh.formCRSMatrix();
				} else {
					A = sch.formCRSMatrix();
				}

				timer.stop("Matrix Formation");

				if (f_m) {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, args::get(f_m).c_str(), FILE_MODE_WRITE,
					                      &viewer);
					MatView(A, viewer);
					PetscViewerDestroy(&viewer);
				}
			}

			if (false) {
#ifdef ENABLE_MUELU
			} else if (f_meulu) {
				timer.start("MeuLu Setup");
				meulusolver = new MueLuWrapper(A, tol, args::get(f_meulu));
				timer.stop("MeuLu Setup");
#endif
#ifdef ENABLE_MUELU_CUDA
			} else if (f_meulucuda) {
				timer.start("MeuLu Setup");
				meulucudasolver = new MueLuCudaWrapper(A, tol, args::get(f_meulucuda));
				timer.stop("MeuLu Setup");
#endif
#ifdef ENABLE_AMGX
			} else if (f_amgx) {
				timer.start("AMGX Setup");
				amgxsolver->setMatrix(A);
				timer.stop("AMGX Setup");
#endif
			} else {
				// preconditoners
				timer.start("Petsc Setup");
				KSPSetOperators(solver, A, A);
				KSPSetUp(solver);
				PC pc;
				KSPGetPC(solver, &pc);
				PCSetUp(pc);
				timer.stop("Petsc Setup");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");
		}
		///////////////////
		// solve start
		///////////////////
		timer.start("Complete Solve");

		if (f_noschur || dc.num_global_domains != 1) {
			timer.start("Linear Solve");
			if (false) {
#ifdef ENABLE_MUELU
			} else if (f_meulu) {
				// solve
				if (f_noschur) {
					meulusolver->solve(u, f);
				} else {
					meulusolver->solve(gamma, b);
				}
#endif
#ifdef ENABLE_MUELU_CUDA
			} else if (f_meulucuda) {
				// solve
				if (f_noschur) {
					meulucudasolver->solve(u, f);
				} else {
					meulucudasolver->solve(gamma, b);
				}
#endif
#ifdef ENABLE_AMGX
			} else if (f_amgx) {
				// solve
				if (f_noschur) {
					amgxsolver->solve(u, f);
				} else {
					amgxsolver->solve(gamma, b);
				}
#endif
			} else {
				KSPSetTolerances(solver, 0.0, tol, PETSC_DEFAULT, 5000);
				if (f_noschur) {
					KSPSolve(solver, f, u);
				} else {
					KSPSolve(solver, b, gamma);
				}
				int its;
				KSPGetIterationNumber(solver, &its);
				if (my_global_rank == 0) {
					cout << "Iterations: " << its << endl;
				}
			}
			timer.stop("Linear Solve");
		}

		if (!f_noschur) {
			// Do one last solve
			timer.start("Patch Solve");

			sch.solveWithInterface(f, u, gamma, diff);

			timer.stop("Patch Solve");
		}

		///////////////////
		// solve end
		///////////////////
		timer.stop("Complete Solve");

		// residual
		PW<Vec> resid = dc.getNewDomainVec();
		PW<Vec> au    = dc.getNewDomainVec();
		if (f_noschur) {
			MatMult(A, u, au);
		} else {
			sch.applyWithInterface(u, gamma, au);
		}
		VecAXPBYPCZ(resid, -1.0, 1.0, 0.0, au, f);
		double residual;
		VecNorm(resid, NORM_2, &residual);
		double fnorm;
		VecNorm(f, NORM_2, &fnorm);

		// error
		PW<Vec> error = dc.getNewDomainVec();
		VecAXPBYPCZ(error, -1.0, 1.0, 0.0, exact, u);
		if (f_neumann) {
			double uavg = dc.integrate(u) / dc.area();
			double eavg = dc.integrate(exact) / dc.area();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			VecShift(error, eavg - uavg);
		}
		double error_norm;
		VecNorm(error, NORM_2, &error_norm);
		double exact_norm;
		VecNorm(exact, NORM_2, &exact_norm);

		double ausum = dc.integrate(au);
		double fsum  = dc.integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = dc.getGlobalNumCells();
			cout << "Total cells: " << total_cells << endl;
		}

		// output
		MMWriter mmwriter(dc, f_amr);
		if (f_s) {
			mmwriter.write(u, args::get(f_s));
		}
		if (f_resid) {
			mmwriter.write(resid, args::get(f_resid));
		}
		if (f_error) {
			mmwriter.write(error, args::get(f_error));
		}
		if (f_g) {
			PetscViewer viewer;
			PetscViewerASCIIOpen(PETSC_COMM_WORLD, args::get(f_g).c_str(), &viewer);
			VecView(gamma, viewer);
		}
		if (f_outclaw) {
			ClawWriter writer(dc);
			writer.write(u, resid);
		}
#ifdef HAVE_VTK
		if (f_outvtk) {
			VtkWriter writer(dc, args::get(f_outvtk));
			writer.add(u, "Solution");
			writer.add(error, "Error");
			writer.add(resid, "Residual");
			writer.write();
		}
#endif
		cout.unsetf(std::ios_base::floatfield);
#ifdef ENABLE_MUELU
		if (meulusolver != nullptr) {
			delete meulusolver;
		}
#endif
#ifdef ENABLE_MUELU_CUDA
		if (meulucudasolver != nullptr) {
			delete meulucudasolver;
		}
#endif
	}

#ifdef ENABLE_AMGX
	if (amgxsolver != nullptr) {
		delete amgxsolver;
	}
#endif
#ifdef ENABLE_MUELU_CUDA
	if (f_meulucuda) {
		MueLuCudaWrapper::finalize();
	}
#endif
	if (my_global_rank == 0) {
		cout << timer;
	}
	PetscFinalize();
	return 0;
}
