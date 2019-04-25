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
#include "FivePtPatchOperator.h"
#include "FunctionWrapper.h"
//#include "IfaceMatrixHelper.h"
#include "BilinearInterpolator.h"
#include "GMG/Helper2d.h"
#include "Init.h"
#include "MatrixHelper2d.h"
#include "PatchSolvers/DftPatchSolver.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "PatchSolvers/FishpackPatchSolver.h"
#include "PolyChebPrec.h"
#include "QuadInterpolator.h"
#include "SchurHelper.h"
#include "SchurMatrixHelper2d.h"
#include "Writers/ClawWriter.h"
#include "Writers/MMWriter.h"
#ifdef HAVE_VTK
#include "Writers/VtkWriter2d.h"
#endif
#ifdef HAVE_P4EST
#include "TreeToP4est.h"
#include "p4estBLG.h"
#endif
#include "CLI11.hpp"
#include "Timer.h"
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
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);

	// parse input
	CLI::App app{"Thunderegg 3d poisson solver example"};

	app.set_config("--config", "", "Read an ini file", false);
	// program options
	int n;
	app.add_option("-n,--num_cells", n, "Number of cells in each direction, on each patch")
	->required();

	int loop_count = 1;
	app.add_option("-l", loop_count, "Number of times to run program");

	string matrix_type = "wrap";
	app.add_set_ignore_case("--matrix_type", matrix_type, {"wrap", "crs"},
	                        "Which type of matrix operator to use");

	int div = 0;
	app.add_option("--divide", div, "Number of levels to add to quadtree");

	bool no_zero_rhs_avg = false;
	app.add_flag("--nozerof", no_zero_rhs_avg,
	             "Make the average of the rhs on neumann problems zero");

	double tolerance = 1e-12;
	app.add_option("-t,--tolerance", tolerance, "Tolerance of Krylov solver");

	bool neumann;
	app.add_flag("--neumann", neumann, "Use neumann boundary conditions");

	bool solve_schur = false;
	app.add_flag("--schur", solve_schur, "Solve the schur compliment system");

	string mesh_filename = "";
	app.add_option("--mesh", mesh_filename, "Filename of mesh to use")
	->required()
	->check(CLI::ExistingFile);

	string problem = "trig";
	app.add_set_ignore_case("--problem", problem, {"trig", "gauss", "zero", "circle"},
	                        "Which problem to solve");

	string preconditioner = "";
	auto   prec_opt
	= app.add_set_ignore_case("--prec", preconditioner, {"GMG"}, "Which Preconditoner to use");

	string gmg_filename = "";
	app.add_option("--gmg_config", gmg_filename, "Which problem to solve")->needs(prec_opt);

	string patch_solver = "fftw";
	app.add_option("--patch_solver", patch_solver, "Which patch solver to use");

	bool setrow = false;
	app.add_flag("--setrow", setrow, "Set row in matrix");

	string petsc_opts = "";
	app.add_option("--petsc_opts", petsc_opts, "petsc options");
	// output options

	string claw_filename = "";
	app.add_option("--out_claw", claw_filename, "Filename of clawpack output");

#ifdef HAVE_VTK
	string vtk_filename = "";
	app.add_option("--out_vtk", vtk_filename, "Filename of vtk output");
#endif

	string matrix_filename = "";
	app.add_option("--out_matrix", matrix_filename, "Filename of matrix output");
	string solution_filename = "";
	app.add_option("--out_solution", solution_filename, "Filename of solution output");
	string resid_filename = "";
	app.add_option("--out_resid", resid_filename, "Filename of residual output");
	string error_filename = "";
	app.add_option("--out_error", error_filename, "Filename of error output");
	string rhs_filename = "";
	app.add_option("--out_rhs", rhs_filename, "Filename of rhs output");
	string gamma_filename = "";
	app.add_option("--out_gamma", gamma_filename, "Filename of gamma output");

#ifdef HAVE_P4EST
	bool use_p4est = false;
	app.add_flag("--p4est", use_p4est, "use p4est");
#endif

	string config_out_filename = "";
	auto   out_config_opt
	= app.add_option("--output_config", config_out_filename, "Save CLI options to config file");

	CLI11_PARSE(app, argc, argv);

	if (config_out_filename != "") {
		app.remove_option(out_config_opt);
		ofstream file_out(config_out_filename);
		file_out << app.config_to_str(true, true);
		file_out.close();
	}

	PetscOptionsInsertString(nullptr, petsc_opts.c_str());

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// Set the number of discretization points in the x and y direction.
	std::array<int, 2> ns;
	ns.fill(n);

	///////////////
	// Create Mesh
	///////////////
	shared_ptr<DomainCollection<2>> dc;
	Tree<2>                         t;
	t = Tree<2>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}
	std::vector<std::map<int, std::shared_ptr<Domain<2>>>> levels;
#ifdef HAVE_P4EST
	if (use_p4est) {
		TreeToP4est ttp(t);
		p4estBLG    blg(ttp.p4est, n);
		levels = blg.levels;
#else
	if (false) {
#endif
	} else {
		BalancedLevelsGenerator<2> blg(t, n);
		// partition domains if running in parallel
		if (num_procs > 1) { blg.zoltanBalance(); }
		levels = blg.levels;
	}

	dc.reset(new DomainCollection<2>(levels[t.num_levels - 1], ns));
	if (neumann) { dc->setNeumann(); }

	// the functions that we are using
	function<double(double, double)> ffun;
	function<double(double, double)> gfun;
	function<double(double, double)> nfun;
	function<double(double, double)> nfuny;

	if (problem == "zero") {
		ffun  = [](double x, double y) { return 0; };
		gfun  = [](double x, double y) { return 0; };
		nfun  = [](double x, double y) { return 0; };
		nfuny = [](double x, double y) { return 0; };
	} else if (problem == "circle") {
		ffun = [](double x, double y) {
			double xdist, ydist, dist;
			// distance form center circle
			xdist = x - 0.5;
			ydist = y - 0.5;
			dist  = sqrt(xdist * xdist + ydist * ydist);
			if (dist < 0.2) { return 1; }
			for (int i = 0; i < 4; i++) {
				// larger circles
				double theta = i * M_PI / 2.0;
				xdist        = x - (0.3 * cos(theta) + 0.5);
				ydist        = y - (0.3 * sin(theta) + 0.5);
				dist         = sqrt(xdist * xdist + ydist * ydist);
				if (dist < 0.1) { return 1; }
				// smaller circles
				theta = M_PI / 4.0 + i * M_PI / 2.0;
				xdist = x - (0.275 * cos(theta) + 0.5);
				ydist = y - (0.275 * sin(theta) + 0.5);
				dist  = sqrt(xdist * xdist + ydist * ydist);
				if (dist < 0.075) { return 1; }
			}
			return 0;
		};
		gfun  = [](double x, double y) { return 0; };
		nfun  = [](double x, double y) { return 0; };
		nfuny = [](double x, double y) { return 0; };
	} else if (problem == "gauss") {
		gfun = [](double x, double y) { return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y)); };
		ffun = [](double x, double y) {
			return 100 * M_PI * M_PI * (pow(sin(10 * M_PI * x), 2) - cos(10 * M_PI * x))
			       * exp(cos(10 * M_PI * x))
			       + 121 * M_PI * M_PI * (cos(11 * M_PI * y) - pow(sin(11 * M_PI * y), 2))
			         * exp(cos(11 * M_PI * y));
		};
		nfun = [](double x, double y) {
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};

		nfuny = [](double x, double y) {
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
	} else {
		ffun
		= [](double x, double y) { return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x); };
		gfun  = [](double x, double y) { return sinl(M_PI * y) * cosl(2 * M_PI * x); };
		nfun  = [](double x, double y) { return -2 * M_PI * sinl(M_PI * y) * sinl(2 * M_PI * x); };
		nfuny = [](double x, double y) { return M_PI * cosl(M_PI * y) * cosl(2 * M_PI * x); };
	}

	// set the patch solver
	shared_ptr<PatchSolver<2>> p_solver;
	if (patch_solver == "dft") {
		p_solver.reset(new DftPatchSolver<2>(*dc));
	} else {
		p_solver.reset(new FftwPatchSolver<2>(*dc));
	}

	// patch operator
	shared_ptr<PatchOperator<2>> p_operator(new FivePtPatchOperator());

	// interface interpolator
	shared_ptr<Interpolator<2>> p_interp(new BilinearInterpolator());

#ifdef ENABLE_AMGX
	AmgxWrapper *amgxsolver = nullptr;
	if (f_amgx) { amgxsolver = new AmgxWrapper(args::get(f_amgx)); }
#endif

#ifdef ENABLE_MUELU_CUDA
	if (f_meulucuda) { MueLuCudaWrapper::initialize(); }
#endif
	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		shared_ptr<SchurHelper<2>> sch(new SchurHelper<2>(dc, p_solver, p_operator, p_interp));
		/*
		if (f_outim) {
		    IfaceMatrixHelper imh(dc);
		    PW<Mat>           IA = imh.formCRSMatrix();
		    PetscViewer       viewer;
		    PetscViewerBinaryOpen(PETSC_COMM_WORLD, args::get(f_outim).c_str(), FILE_MODE_WRITE,
		                          &viewer);
		    MatView(IA, viewer);
		    PetscViewerDestroy(&viewer);
		}
		*/

		shared_ptr<PetscVector<2>> u     = dc->getNewDomainVec();
		shared_ptr<PetscVector<2>> exact = dc->getNewDomainVec();
		shared_ptr<PetscVector<2>> f     = dc->getNewDomainVec();

		if (neumann) {
			Init::initNeumann2d(*dc, n, f->vec, exact->vec, ffun, gfun, nfun, nfuny);
		} else {
			Init::initDirichlet2d(*dc, n, f->vec, exact->vec, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		// Create the gamma and diff vectors
		shared_ptr<PetscVector<1>> gamma = sch->getNewSchurVec();
		shared_ptr<PetscVector<1>> diff  = sch->getNewSchurVec();
		shared_ptr<PetscVector<1>> b     = sch->getNewSchurVec();
		PW<Mat>                    A;
		shared_ptr<GMG::Helper2d>  gh;

		// Create linear problem for the Belos solver
		PW<KSP> solver;
		KSPCreate(MPI_COMM_WORLD, &solver);
		KSPSetFromOptions(solver);

		if (neumann && !no_zero_rhs_avg) {
			double fdiff = dc->integrate(f) / dc->volume();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
			f->shift(-fdiff);
		}

#ifdef ENABLE_MUELU
		MueLuWrapper *meulusolver = nullptr;
#endif
#ifdef ENABLE_MUELU_CUDA
		MueLuCudaWrapper *meulucudasolver = nullptr;
#endif

		if (!solve_schur || dc->num_global_domains != 1) {
			// do iterative solve

			if (solve_schur) {
				// Get the b vector
				gamma->set(0);
				sch->solveWithInterface(f, u, gamma, b);
				b->scale(-1);

				if (rhs_filename != "") {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, rhs_filename.c_str(), FILE_MODE_WRITE,
					                      &viewer);
					VecView(b->vec, viewer);
					PetscViewerDestroy(&viewer);
				}
			}

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			timer.start("Matrix Formation");
			if (matrix_type == "wrap") {
				// Create a function wrapper
				if (solve_schur) {
					A = FuncWrap<2>::getMatrix(sch.get(), dc.get());
				} else {
					A = FullFuncWrap<2>::getMatrix(sch.get(), dc.get());
				}
			} else {
				if (solve_schur) {
					SchurMatrixHelper2d mh(sch);
					A = mh.formCRSMatrix();
				} else {
					MatrixHelper2d mh(*dc);
					A = mh.formCRSMatrix();
				}
				if (setrow) {
					int row = 0;
					MatZeroRows(A, 1, &row, 1.0, nullptr, nullptr);
				}
			}
			timer.stop("Matrix Formation");

			if (matrix_filename != "") {
				PW<Mat> A_copy;
				MatConvert(A, MATMPIAIJ, MAT_INITIAL_MATRIX, &A_copy);
				PetscViewer viewer;
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_filename.c_str(), FILE_MODE_WRITE,
				                      &viewer);
				MatView(A_copy, viewer);
				PetscViewerDestroy(&viewer);
			}

			// preconditoners
			timer.start("Petsc Setup");
			KSPSetOperators(solver, A, A);
			KSPSetUp(solver);
			PC pc;
			KSPGetPC(solver, &pc);
			if (preconditioner == "GMG") {
				vector<shared_ptr<DomainCollection<2>>> dcs(t.num_levels);
				dcs[0] = dc;
				for (int i = 1; i < t.num_levels; i++) {
					dcs[i].reset(new DomainCollection<2>(levels[t.num_levels - 1 - i], ns));
				}
				gh.reset(new GMG::Helper2d(n, dcs, sch, gmg_filename));
				gh->getPrec(pc);
			}
			PCSetUp(pc);
			timer.stop("Petsc Setup");
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");
		}
		///////////////////
		// solve start
		///////////////////
		timer.start("Complete Solve");

		if (!solve_schur || dc->num_global_domains != 1) {
			timer.start("Linear Solve");
			KSPSetTolerances(solver, tolerance, tolerance, PETSC_DEFAULT, 5000);
			if (solve_schur) {
				KSPSolve(solver, b->vec, gamma->vec);
			} else {
				KSPSolve(solver, f->vec, u->vec);
			}
			int its;
			KSPGetIterationNumber(solver, &its);
			if (my_global_rank == 0) { cout << "Iterations: " << its << endl; }
			timer.stop("Linear Solve");
		}

		if (solve_schur) {
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
		shared_ptr<PetscVector<2>> resid = dc->getNewDomainVec();
		shared_ptr<PetscVector<2>> au    = dc->getNewDomainVec();
		if (solve_schur) {
			sch->applyWithInterface(u, gamma, au);
		} else {
			MatMult(A, u->vec, au->vec);
		}
		VecAXPBYPCZ(resid->vec, -1.0, 1.0, 0.0, au->vec, f->vec);
		double residual = resid->twoNorm();
		double fnorm    = f->twoNorm();

		// error
		shared_ptr<PetscVector<2>> error = dc->getNewDomainVec();
		VecAXPBYPCZ(error->vec, -1.0, 1.0, 0.0, exact->vec, u->vec);
		if (neumann) {
			double uavg = dc->integrate(u) / dc->volume();
			double eavg = dc->integrate(exact) / dc->volume();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			VecShift(error->vec, eavg - uavg);
		}
		double error_norm = error->twoNorm();
		double exact_norm = exact->twoNorm();

		double ausum = dc->integrate(au);
		double fsum  = dc->integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = dc->getGlobalNumCells();
			cout << "Total cells: " << total_cells << endl;
		}

		// output
		if (claw_filename != "") {
			ClawWriter writer(*dc);
			writer.write(u->vec, resid->vec);
		}
#ifdef HAVE_VTK
		if (vtk_filename != "") {
			VtkWriter2d writer(*dc, vtk_filename);
			writer.add(u->vec, "Solution");
			writer.add(error->vec, "Error");
			writer.add(resid->vec, "Residual");
			writer.add(f->vec, "RHS");
			writer.add(exact->vec, "Exact");
			writer.write();
		}
#endif
		cout.unsetf(std::ios_base::floatfield);
#ifdef ENABLE_MUELU
		if (meulusolver != nullptr) { delete meulusolver; }
#endif
#ifdef ENABLE_MUELU_CUDA
		if (meulucudasolver != nullptr) { delete meulucudasolver; }
#endif
	}

#ifdef ENABLE_AMGX
	if (amgxsolver != nullptr) { delete amgxsolver; }
#endif
#ifdef ENABLE_MUELU_CUDA
	if (f_meulucuda) { MueLuCudaWrapper::finalize(); }
#endif
	if (my_global_rank == 0) { cout << timer; }
	PetscFinalize();
	return 0;
}
