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

#include "Init.h"
#include "Writers/ClawWriter.h"
#include <Thunderegg/BiCGStab.h>
#include <Thunderegg/BilinearInterpolator.h>
#include <Thunderegg/Domain.h>
#include <Thunderegg/FivePtPatchOperator.h>
#include <Thunderegg/GMG/CycleFactory2d.h>
#include <Thunderegg/MatrixHelper2d.h>
#include <Thunderegg/Operators/DomainWrapOp.h>
#include <Thunderegg/Operators/PetscMatOp.h>
#include <Thunderegg/Operators/SchurWrapOp.h>
#include <Thunderegg/PatchSolvers/DftPatchSolver.h>
#include <Thunderegg/PatchSolvers/FftwPatchSolver.h>
#include <Thunderegg/PetscShellCreator.h>
#include <Thunderegg/PolyChebPrec.h>
#include <Thunderegg/QuadInterpolator.h>
#include <Thunderegg/SchurHelper.h>
#include <Thunderegg/SchurMatrixHelper2d.h>
#include <Thunderegg/SchwarzPrec.h>
#include <Thunderegg/StarPatchOp.h>
#include <Thunderegg/ThundereggDomGen.h>
#include <Thunderegg/Timer.h>
#ifdef HAVE_VTK
#include "Writers/VtkWriter2d.h"
#endif
#ifdef HAVE_P4EST
#include "TreeToP4est.h"
#include <Thunderegg/P4estDomGen.h>
#endif
#include "CLI11.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <string>
#include <unistd.h>
//#include "IfaceMatrixHelper.h"

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

	string solver_type = "thunderegg";
	app.add_set_ignore_case("--solver", solver_type, {"petsc", "thunderegg"},
	                        "Which Solver to use");

	string preconditioner = "";
	app.add_set_ignore_case("--prec", preconditioner, {"GMG"}, "Which Preconditoner to use");

	string patch_solver = "fftw";
	app.add_option("--patch_solver", patch_solver, "Which patch solver to use");

	bool setrow = false;
	app.add_flag("--setrow", setrow, "Set row in matrix");

	string petsc_opts = "";
	app.add_option("--petsc_opts", petsc_opts, "petsc options");

	// GMG options

	auto gmg = app.add_subcommand("GMG", "GMG solver options");

	GMG::CycleOpts copts;

	gmg->add_option("--max_levels", copts.max_levels,
	                "The max number of levels in GMG cycle. 0 means no limit.");

	gmg->add_option(
	"--patches_per_proc", copts.patches_per_proc,
	"Lowest level is guaranteed to have at least this number of patches per processor.");

	gmg->add_option("--pre_sweeps", copts.pre_sweeps, "Number of sweeps on down cycle");

	gmg->add_option("--post_sweeps", copts.post_sweeps, "Number of sweeps on up cycle");

	gmg->add_option("--mid_sweeps", copts.mid_sweeps,
	                "Number of sweeps inbetween up and down cycle");

	gmg->add_option("--coarse_sweeps", copts.coarse_sweeps, "Number of sweeps on coarse level");

	gmg->add_option("--cycle_type", copts.cycle_type, "Cycle type");

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
	shared_ptr<Domain<2>> domain;
	Tree<2>               t;
	t = Tree<2>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}

	shared_ptr<DomainGenerator<2>> dcg;
#ifdef HAVE_P4EST
	if (use_p4est) {
		TreeToP4est ttp(t);

		auto bmf = [](int block_no, double unit_x, double unit_y, double &x, double &y) {
			x = unit_x;
			y = unit_y;
		};
		auto inf
		= [=](Side<2> s, const array<double, 2> &, const array<double, 2> &) { return neumann; };

		dcg.reset(new P4estDomGen(ttp.p4est, ns, inf, bmf));
#else
	if (false) {
#endif
	} else {
		dcg.reset(new ThundereggDomGen<2>(t, ns, neumann));
	}

	domain = dcg->getFinestDomain();

	// the functions that we are using
	function<double(double, double)> ffun;
	function<double(double, double)> gfun;
	function<double(double, double)> nfun;
	function<double(double, double)> nfuny;

	if (problem == "gauss") {
		double x0    = .5;
		double y0    = .5;
		double alpha = 1000;
		ffun         = [&](double x, double y) {
            double r2 = pow((x - x0), 2) + pow((y - y0), 2);
            return exp(-alpha / 2.0 * r2) * (pow(alpha, 2) * r2 - 2 * alpha);
		};
		gfun = [&](double x, double y) {
			double r2 = pow((x - x0), 2) + pow((y - y0), 2);
			return exp(-alpha / 2.0 * r2);
		};
		nfun  = [](double x, double y) { return 0; };
		nfuny = [](double x, double y) { return 0; };
	} else if (problem == "zero") {
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
	} else if (problem == "trig gauss") {
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
		p_solver.reset(new DftPatchSolver<2>(*domain));
	} else {
		p_solver.reset(new FftwPatchSolver<2>(*domain));
	}

	// patch operator
	shared_ptr<PatchOperator<2>> p_operator(new StarPatchOp<2>());

	// interface interpolator
	shared_ptr<IfaceInterp<2>> p_interp(new BilinearInterpolator());

	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		shared_ptr<SchurHelper<2>> sch(new SchurHelper<2>(domain, p_solver, p_operator, p_interp));
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

		shared_ptr<PetscVector<2>> u     = domain->getNewDomainVec();
		shared_ptr<PetscVector<2>> exact = domain->getNewDomainVec();
		shared_ptr<PetscVector<2>> f     = domain->getNewDomainVec();
		shared_ptr<PetscVector<2>> au    = domain->getNewDomainVec();

		if (neumann) {
			Init::initNeumann2d(*domain, f->vec, exact->vec, ffun, gfun, nfun, nfuny);
		} else {
			Init::initDirichlet2d(*domain, f->vec, exact->vec, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		// Create the gamma and diff vectors
		shared_ptr<PetscVector<1>> gamma = sch->getNewSchurVec();
		shared_ptr<PetscVector<1>> diff  = sch->getNewSchurVec();
		shared_ptr<PetscVector<1>> b     = sch->getNewSchurVec();

		// Create linear problem for the Belos solver
		PW<KSP> solver;
		KSPCreate(MPI_COMM_WORLD, &solver);
		KSPSetFromOptions(solver);

		if (neumann && !no_zero_rhs_avg) {
			double fdiff = domain->integrate(f) / domain->volume();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
			f->shift(-fdiff);
		}

		if (solve_schur) {
			// initialize vectors for Schur compliment system

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
			std::shared_ptr<Operator<1>> A;
			PW<Mat>                      A_petsc;
			std::shared_ptr<Operator<1>> M;
			PW<PC>                       M_petsc;

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			timer.start("Matrix Formation");
			if (matrix_type == "wrap") {
				A.reset(new SchurWrapOp<2>(domain, sch));
			} else if (matrix_type == "crs") {
				SchurMatrixHelper2d smh(sch);
				A_petsc = smh.formCRSMatrix();
				A.reset(new PetscMatOp<1>(A_petsc));
				if (setrow) {
					int row = 0;
					MatZeroRows(A_petsc, 1, &row, 1.0, nullptr, nullptr);
				}

				if (matrix_filename != "") {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_filename.c_str(),
					                      FILE_MODE_WRITE, &viewer);
					MatView(A_petsc, viewer);
					PetscViewerDestroy(&viewer);
				}
			}
			timer.stop("Matrix Formation");
			// preconditoners
			timer.start("Preconditioner Setup");
			if (preconditioner == "Scwharz") {
				throw 3;
			} else if (preconditioner == "GMG") {
				throw 3;
			}
			timer.stop("Preconditioner Setup");

			PW<KSP> solver;
			// setup petsc if needed
			if (solver_type == "petsc") {
				timer.start("Petsc Setup");

				KSPCreate(MPI_COMM_WORLD, &solver);
				KSPSetFromOptions(solver);
				KSPSetOperators(solver, A_petsc, A_petsc);
				if (M != nullptr) {
					PC M_petsc;
					KSPGetPC(solver, &M_petsc);
					PetscShellCreator::getPCShell(M_petsc, M, sch);
				}
				KSPSetUp(solver);
				KSPSetTolerances(solver, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, 5000);

				timer.stop("Petsc Setup");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");

			timer.start("Linear Solve");
			if (solver_type == "petsc") {
				KSPSolve(solver, b->vec, gamma->vec);
				int its;
				KSPGetIterationNumber(solver, &its);
				if (my_global_rank == 0) { cout << "Iterations: " << its << endl; }
			} else {
				std::shared_ptr<VectorGenerator<1>> vg(new SchurHelperVG<1>(sch));

				int its = BiCGStab<1>::solve(vg, A, gamma, b, M);
				if (my_global_rank == 0) { cout << "Iterations: " << its << endl; }
			}
			timer.stop("Linear Solve");

			// Do one last solve
			timer.start("Patch Solve");

			sch->solveWithInterface(f, u, gamma, diff);

			timer.stop("Patch Solve");

			sch->applyWithInterface(u, gamma, au);
		} else {
			std::shared_ptr<Operator<2>> A;
			PW<Mat>                      A_petsc;
			std::shared_ptr<Operator<2>> M;
			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			timer.start("Matrix Formation");
			if (matrix_type == "wrap") {
				A.reset(new DomainWrapOp<2>(sch));
				A_petsc = PetscShellCreator::getMatShell(A, domain);
			} else if (matrix_type == "crs") {
				MatrixHelper2d mh(domain);
				A_petsc = mh.formCRSMatrix();
				A.reset(new PetscMatOp<2>(A_petsc));
				if (setrow) {
					int row = 0;
					MatZeroRows(A_petsc, 1, &row, 1.0, nullptr, nullptr);
				}

				if (matrix_filename != "") {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_filename.c_str(),
					                      FILE_MODE_WRITE, &viewer);
					MatView(A_petsc, viewer);
					PetscViewerDestroy(&viewer);
				}
			} else if (matrix_type == "pbm") {
				throw 3;
			}
			timer.stop("Matrix Formation");
			// preconditoners
			timer.start("Preconditioner Setup");
			if (preconditioner == "Scwharz") {
				M.reset(new SchwarzPrec<2>(sch));
			} else if (preconditioner == "GMG") {
				timer.start("GMG Setup");

				M = GMG::CycleFactory2d::getCycle(copts, dcg, p_solver, p_operator, p_interp);

				timer.stop("GMG Setup");
			} else if (preconditioner == "cheb") {
				throw 3;
			}
			timer.stop("Preconditioner Setup");

			PW<KSP> solver;
			// setup petsc if needed
			if (solver_type == "petsc") {
				timer.start("Petsc Setup");

				KSPCreate(MPI_COMM_WORLD, &solver);
				KSPSetFromOptions(solver);
				KSPSetOperators(solver, A_petsc, A_petsc);
				if (M != nullptr) {
					PC M_petsc;
					KSPGetPC(solver, &M_petsc);
					PetscShellCreator::getPCShell(M_petsc, M, domain);
				}
				KSPSetUp(solver);
				KSPSetTolerances(solver, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, 5000);

				timer.stop("Petsc Setup");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");

			timer.start("Linear Solve");
			if (solver_type == "petsc") {
				KSPSolve(solver, f->vec, u->vec);
				int its;
				KSPGetIterationNumber(solver, &its);
				if (my_global_rank == 0) { cout << "Iterations: " << its << endl; }
			} else {
				std::shared_ptr<VectorGenerator<2>> vg(new DomainVG<2>(domain));

				int its = BiCGStab<2>::solve(vg, A, u, f, M);
				if (my_global_rank == 0) { cout << "Iterations: " << its << endl; }
			}
			timer.stop("Linear Solve");

			A->apply(u, au);
		}

		// residual
		shared_ptr<PetscVector<2>> resid = domain->getNewDomainVec();
		VecAXPBYPCZ(resid->vec, -1.0, 1.0, 0.0, au->vec, f->vec);
		double residual = resid->twoNorm();
		double fnorm    = f->twoNorm();

		// error
		shared_ptr<PetscVector<2>> error = domain->getNewDomainVec();
		VecAXPBYPCZ(error->vec, -1.0, 1.0, 0.0, exact->vec, u->vec);
		if (neumann) {
			double uavg = domain->integrate(u) / domain->volume();
			double eavg = domain->integrate(exact) / domain->volume();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			VecShift(error->vec, eavg - uavg);
		}
		double error_norm = error->twoNorm();
		double exact_norm = exact->twoNorm();

		double ausum = domain->integrate(au);
		double fsum  = domain->integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = domain->getNumGlobalCells();
			cout << "Total cells: " << total_cells << endl;
		}

		// output
		if (claw_filename != "") {
			ClawWriter writer(domain);
			writer.write(u->vec, resid->vec);
		}
#ifdef HAVE_VTK
		if (vtk_filename != "") {
			VtkWriter2d writer(*domain, vtk_filename);
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
