/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively 
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#include "BalancedLevelsGenerator.h"
#include "DomainCollection.h"
#include "FunctionWrapper.h"
#include "GMG/Helper.h"
#include "Init.h"
#include "MatrixHelper.h"
#include "OctTree.h"
#include "PatchSolvers/DftPatchSolver.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "PolyChebPrec.h"
#include "SchurHelper.h"
#include "SchurMatrixHelper.h"
#include "SevenPtPatchOperator.h"
#include "Timer.h"
#include "TriLinInterp.h"
#include "Writers/ClawWriter.h"
#include "Writers/MMWriter.h"
#include "args.hxx"
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
	args::ValueFlag<int>    f_cube(parser, "num_domains", "create a num_domains^3 cube of grids",
                                {"cube"});
	args::ValueFlag<int>    f_amr(parser, "num_domains",
                               "create a num_domains x num_domains square "
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
	args::Flag              f_scharz(parser, "", "use schwarz preconditioner", {"schwarz"});
	args::ValueFlag<string> f_gmg(parser, "config_file", "use GMG preconditioner", {"gmg"});
	args::Flag              f_cfft(parser, "", "use GMG preconditioner", {"cfft"});
	args::Flag              f_pbm(parser, "", "use GMG preconditioner", {"pbm"});
	args::Flag              f_ibd(parser, "", "use GMG preconditioner", {"ibd"});
	args::Flag              f_cheb(parser, "", "cheb preconditioner", {"cheb"});
	args::Flag              f_dft(parser, "", "dft", {"dft"});

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
	if (f_t) { tol = args::get(f_t); }

	int loop_count = 1;
	if (f_l) { loop_count = args::get(f_l); }
	/***********
	 * Input parsing done
	 **********/

	///////////////
	// Create Mesh
	///////////////
	shared_ptr<DomainCollection<3>> dc;
	Tree<3>                         t;
	if (f_mesh) {
		string d = args::get(f_mesh);
		t        = Tree<3>(d);
	}
	if (f_div) {
		for (int i = 0; i < args::get(f_div); i++) {
			t.refineLeaves();
		}
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
			return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y)) + exp(cos(12 * M_PI * z));
		};
		ffun = [](double x, double y, double z) {
            constexpr double pi2 = M_PI*M_PI;
			return
                -M_PI * M_PI
			       * (100 * exp(cos(10 * M_PI * x)) * cos(10 * M_PI * x)
			          - 100 * exp(cos(10 * M_PI * x)) * pow(sin(10 * M_PI * x), 2)
			          - 121 * exp(cos(11 * M_PI * y)) * cos(11 * M_PI * y)
			          + 121 * exp(cos(11 * M_PI * y)) * pow(sin(11 * M_PI * y), 2)
			          + 144 * exp(cos(12 * M_PI * z)) * cos(12 * M_PI * z)
			          - 144 * exp(cos(12 * M_PI * z)) * pow(sin(12 * M_PI * z), 2));
		};
		nfunx = [](double x, double y, double z) {
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};
		nfuny = [](double x, double y, double z) {
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
		nfunz = [](double x, double y, double z) {
			return -12 * M_PI * sin(12 * M_PI * z) * exp(cos(12 * M_PI * z));
		};
	} else {
		ffun = [](double x, double y, double z) {
            x+=.3;
            y+=.3;
            z+=.3;
			return -77.0 / 36 * M_PI * M_PI * sin(M_PI * x) * cos(2.0 / 3 * M_PI * y)
			       * sin(5.0 / 6 * M_PI * z);
		};
		gfun = [](double x, double y, double z) {
            x+=.3;
            y+=.3;
            z+=.3;
			return sin(M_PI * x) * cos(2.0 / 3 * M_PI * y) * sin(5.0 / 6 * M_PI * z);
		};
		nfunx = [](double x, double y, double z) {
            x+=.3;
            y+=.3;
            z+=.3;
			return M_PI * cos(M_PI * x) * cos(2.0 / 3 * M_PI * y) * sin(5.0 / 6 * M_PI * z);
		};
		nfuny = [](double x, double y, double z) {
            x+=.3;
            y+=.3;
            z+=.3;
			return -2.0 / 3 * M_PI * sin(M_PI * x) * sin(2.0 / 3 * M_PI * y)
			       * sin(5.0 / 6 * M_PI * z);
		};
		nfunz = [](double x, double y, double z) {
            x+=.3;
            y+=.3;
            z+=.3;
			return 5.0 / 6 * M_PI * sin(M_PI * x) * cos(2.0 / 3 * M_PI * y)
			       * cos(5.0 / 6 * M_PI * z);
		};
	}

	// set the patch solver
	shared_ptr<PatchSolver<3>> p_solver;

	// patch operator
	shared_ptr<PatchOperator<3>> p_operator(new SevenPtPatchOperator());

	// interface interpolator
	shared_ptr<Interpolator<3>> p_interp(new TriLinInterp());

#ifdef ENABLE_AMGX
	AmgxWrapper *amgxsolver = nullptr;
	if (f_amgx) { amgxsolver = new AmgxWrapper(args::get(f_amgx)); }
#endif

	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");
		BalancedLevelsGenerator<3> blg(t, n);

		// partition domains if running in parallel
		if (num_procs > 1) {
			timer.start("Zoltan Balance");
			blg.zoltanBalance();
			timer.stop("Zoltan Balance");
		}

		dc.reset(new DomainCollection<3>(blg.levels[t.num_levels - 1], n));
		if (f_neumann) { dc->setNeumann(); }

		if (f_dft) {
			p_solver.reset(new DftPatchSolver<3>(*dc));
		} else {
			p_solver.reset(new FftwPatchSolver<3>(*dc));
		}
		shared_ptr<SchurHelper<3>> sch(new SchurHelper<3>(*dc, p_solver, p_operator, p_interp));
		MatrixHelper               mh(*dc);

		PW<Vec> u     = dc->getNewDomainVec();
		PW<Vec> exact = dc->getNewDomainVec();
		PW<Vec> f     = dc->getNewDomainVec();

		if (f_neumann) {
			Init::initNeumann(*dc, n, f, exact, ffun, gfun, nfunx, nfuny, nfunz);
		} else {
			Init::initDirichlet(*dc, n, f, exact, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		// Create the gamma and diff vectors
		PW<Vec>                 gamma = sch->getNewSchurVec();
		PW<Vec>                 diff  = sch->getNewSchurVec();
		PW<Vec>                 b     = sch->getNewSchurVec();
		PW<Mat>                 A;
		shared_ptr<FuncWrap<3>> w;
		shared_ptr<SchwarzPrec> sp;
		shared_ptr<GMG::Helper> gh;

		// Create linear problem for the Belos solver
		PW<KSP> solver;
		KSPCreate(MPI_COMM_WORLD, &solver);
		KSPSetFromOptions(solver);

		if (f_neumann && !f_nozerof) {
			double fdiff = dc->integrate(f) / dc->volume();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
			VecShift(f, -fdiff);
		}

		if (f_noschur || dc->num_global_domains != 1) {
			// do iterative solve

			if (!f_noschur) {
				// Get the b vector
				VecScale(gamma, 0);
				sch->solveWithInterface(f, u, gamma, b);
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
				if (f_noschur) {
					A = FullFuncWrap<3>::getMatrix(sch.get(), dc.get());
				} else {
					A = FuncWrap<3>::getMatrix(sch.get(), dc.get());
				}
			} else {
				timer.start("Matrix Formation");

				if (f_noschur) {
					A = mh.formCRSMatrix();
				} else {
					SchurMatrixHelper smh(sch);
					if (f_pbm) {
						A = smh.getPBMatrix();
					} else {
						A = smh.formCRSMatrix();
					}
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
				if (f_scharz) {
					sp.reset(new SchwarzPrec(sch.get(), &*dc));
					sp->getPrec(pc);
				}
				if (f_gmg) {
					timer.start("GMG Setup");
					timer.start("GMG Domain Collection Setup");
					vector<shared_ptr<DomainCollection<3>>> dcs(t.num_levels);
					dcs[0] = dc;
					for (int i = 1; i < t.num_levels; i++) {
						dcs[i].reset(new DomainCollection<3>(blg.levels[t.num_levels - 1 - i], n));
					}
					timer.stop("GMG Domain Collection Setup");

					gh.reset(new GMG::Helper(n, dcs, sch, args::get(f_gmg)));
					timer.stop("GMG Setup");
					gh->getPrec(pc);
				}
				if (f_cheb) {
					PolyChebPrec *pcp = new PolyChebPrec(*sch, *dc);
					pcp->getPrec(pc);
					PCSetUp(pc);
				}
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

		if ((f_noschur || dc->num_global_domains != 1) && !f_cfft) {
			timer.start("Linear Solve");
			if (false) {
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
				KSPSetTolerances(solver, tol, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
				if (f_noschur) {
					KSPSolve(solver, f, u);
				} else {
					KSPSolve(solver, b, gamma);
				}
				int its;
				KSPGetIterationNumber(solver, &its);
				if (my_global_rank == 0) { cout << "Iterations: " << its << endl; }
			}
			timer.stop("Linear Solve");
		}

		if (!f_noschur) {
			// Do one last solve
			timer.start("Patch Solve");

			sch->solveWithInterface(f, u, gamma, diff);

			timer.stop("Patch Solve");
		}

		///////////////////
		// solve end
		///////////////////
		timer.stop("Complete Solve");

		// residual
		PW<Vec> resid = dc->getNewDomainVec();
		PW<Vec> au    = dc->getNewDomainVec();
		if (f_noschur) {
			MatMult(A, u, au);
		} else {
			sch->applyWithInterface(u, gamma, au);
		}
		VecAXPBYPCZ(resid, -1.0, 1.0, 0.0, au, f);
		double residual;
		VecNorm(resid, NORM_2, &residual);
		double fnorm;
		VecNorm(f, NORM_2, &fnorm);

		// error
		PW<Vec> error = dc->getNewDomainVec();
		VecAXPBYPCZ(error, -1.0, 1.0, 0.0, exact, u);
		if (f_neumann) {
			double uavg = dc->integrate(u) / dc->volume();
			double eavg = dc->integrate(exact) / dc->volume();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			VecShift(error, eavg - uavg);
		}
		double error_norm;
		VecNorm(error, NORM_2, &error_norm);
		double error_norm_inf;
		VecNorm(error, NORM_INFINITY, &error_norm_inf);
		double exact_norm;
		VecNorm(exact, NORM_2, &exact_norm);

		double ausum = dc->integrate(au);
		double fsum  = dc->integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error (2-norm):   " << error_norm / exact_norm << endl;
			std::cout << "Error (inf-norm): " << error_norm_inf << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = dc->getGlobalNumCells();
			cout << "Total cells: " << total_cells << endl;
			cout << "Cores: " << num_procs << endl;
		}

		// output
		if (f_g) {
			PetscViewer viewer;
			PetscViewerASCIIOpen(PETSC_COMM_WORLD, args::get(f_g).c_str(), &viewer);
			VecView(gamma, viewer);
		}
#ifdef HAVE_VTK
		if (f_outvtk) {
			VtkWriter writer(*dc, args::get(f_outvtk));
			writer.add(u, "Solution");
			writer.add(error, "Error");
			writer.add(exact, "Exact");
			writer.add(resid, "Residual");
			writer.add(f, "RHS");
			writer.write();
		}
#endif
		cout.unsetf(std::ios_base::floatfield);
	}

#ifdef ENABLE_AMGX
	if (amgxsolver != nullptr) { delete amgxsolver; }
#endif
	if (my_global_rank == 0) { cout << timer; }
	PetscFinalize();
	return 0;
}
