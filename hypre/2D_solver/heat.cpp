#include "DomainCollection.h"
#include "DomainSignatureCollection.h"
#include "args.h"
#include <HYPRE_krylov.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <valarray>
//#include <mpi.h>
#include <string>
#include <unistd.h>
#include <fstream>


// =========== //
// main driver //
// =========== //

using namespace std;
int main(int argc, char *argv[])
{
	//    using Belos::FuncWrap;
	using namespace std::chrono;

	MPI_Init(NULL, NULL);
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

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
	args::ValueFlag<string> f_read_gamma(parser, "gamma filename",
	                                     "the file to read gamma vector from", {"readgamma"});
	args::ValueFlag<string> f_flux(parser, "flux filename", "the file to write flux difference to",
	                            {"flux"});
	args::ValueFlag<string> f_p(parser, "preconditioner filename",
	                            "the file to write the preconditioner to", {'p'});
	args::ValueFlag<double> f_t(
	parser, "tolerance", "set the tolerance of the iterative solver (default is 1e-10)", {'t'});
    args::ValueFlag<int> f_d(
	parser, "row", "pin gamma value to zero (by modifying that row of the schur compliment matrix)",
	{'z'});
	args::Flag f_pinv(parser, "wrapper", "compute using pseudoinverse", {"pinv"});
	args::Flag f_wrapper(parser, "wrapper", "use a function wrapper", {"wrap"});
	args::Flag f_gauss(parser, "gauss", "solve gaussian function", {"gauss"});
	args::Flag f_prec(parser, "prec", "use block diagonal preconditioner", {"prec"});
	args::Flag f_precata(parser, "prec", "use block diagonal preconditioner", {"precata"});
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {"neumann"});
	args::Flag f_cg(parser, "gmres", "use CG for iterative solver", {"cg"});
	args::Flag f_gmres(parser, "gmres", "use BiCGSTAB for iterative solver", {"gmres"});
	args::Flag f_lsqr(parser, "gmres", "use BiCGSTAB for iterative solver", {"lsqr"});
	args::Flag f_rgmres(parser, "rgmres", "use GCRO-DR (Recycling BiCGSTAB) for iterative solver",
	                    {"rgmres"});
	args::Flag f_bicg(parser, "gmres", "use BiCGStab for iterative solver", {"bicg"});
	args::Flag f_nozero(parser, "nozero", "don't make the average of vector zero in CG solver",
	                    {"nozero"});
	args::Flag f_nozerou(parser, "zerou", "don't modify make so that it zeros the solution",
	                    {"nozerou"});
	args::Flag f_nozerof(parser, "zerou", "don't modify make so that it zeros the solution",
	                    {"nozerof"});
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
		dsc      = DomainSignatureCollection(d, my_global_rank);
	} else if (f_amr) {
		int d = args::get(f_amr);
		dsc   = DomainSignatureCollection(d, d, my_global_rank, true);
	} else {
        int d = args::get(f_square);
		dsc = DomainSignatureCollection(d, d, my_global_rank);
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
	double tol = 1e-10;
	if (f_t) {
		tol = args::get(f_t);
	}

	int loop_count = 1;
	if (f_l) {
		loop_count = args::get(f_l);
	}

    int del = -1;
	if (f_d) {
		del = args::get(f_d);
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
		gfun
		= [](double x, double y) { return exp(cos(10 * M_PIl * x)) - exp(cos(11 * M_PIl * y)); };
		ffun = [](double x, double y) {
			return 100 * M_PIl * M_PIl * (pow(sin(10 * M_PIl * x), 2) - cos(10 * M_PIl * x))
			       * exp(cos(10 * M_PIl * x))
			       + 121 * M_PIl * M_PIl * (cos(11 * M_PIl * y) - pow(sin(11 * M_PIl * y), 2))
			         * exp(cos(11 * M_PIl * y));
		};
		nfunx = [](double x, double y) {
			return -10 * M_PIl * sin(10 * M_PIl * x) * exp(cos(10 * M_PIl * x));
		};

		nfuny = [](double x, double y) {
			return 11 * M_PIl * sin(11 * M_PIl * y) * exp(cos(11 * M_PIl * y));
		};
	} else {
		ffun = [](double x, double y) {
			return -5 * M_PIl * M_PIl * sinl(M_PIl * y) * cosl(2 * M_PIl * x);
		};
		gfun = [](double x, double y) { return sinl(M_PIl * y) * cosl(2 * M_PIl * x); };
		nfunx
		= [](double x, double y) { return -2 * M_PIl * sinl(M_PIl * y) * sinl(2 * M_PIl * x); };
		nfuny = [](double x, double y) { return M_PIl * cosl(M_PIl * y) * cosl(2 * M_PIl * x); };
	}

	valarray<double> times(loop_count);
	for (int loop = 0; loop < loop_count; loop++) {
        MPI_Barrier(MPI_COMM_WORLD);

		steady_clock::time_point domain_start = steady_clock::now();

		DomainCollection dc(dsc, nx);

		if (f_neumann) {
			dc.initNeumann(ffun, gfun, nfunx, nfuny, f_amr);
		} else {
			dc.initDirichlet(ffun, gfun);
            dc.amr = f_amr;
		}


        MPI_Barrier(MPI_COMM_WORLD);
		steady_clock::time_point domain_stop = steady_clock::now();
		duration<double>         domain_time = domain_stop - domain_start;

		if (my_global_rank == 0)
			cout << "Domain Initialization Time: " << domain_time.count() << endl;

		if (f_neumann && !f_nozerof) {
			double fdiff = (dc.integrateBoundaryFlux() - dc.integrateF()) / dc.area();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
			for (auto &p : dc.domains) {
				Domain &d = p.second;
				d.f += fdiff;
			}
		}

		//************
		// SOLVE
		//************
		// initialize the x and b vectors
		dc.initVectors();
		HYPRE_SStructSolver solver;

		HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
		HYPRE_SStructBiCGSTABSetMaxIter(solver, 5000);
		HYPRE_SStructBiCGSTABSetTol(solver, tol);
		HYPRE_SStructBiCGSTABSetPrintLevel(solver, 3);
		HYPRE_SStructBiCGSTABSetLogging(solver, 1);

		// A, x, and b are stored in the DomainCollection object
		HYPRE_SStructBiCGSTABSetup(solver, dc.A, dc.b, dc.x);
		HYPRE_SStructBiCGSTABSolve(solver, dc.A, dc.b, dc.x);

		int    num_iterations;
		double final_res_norm;
		HYPRE_SStructBiCGSTABGetNumIterations(solver, &num_iterations);
		HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

		HYPRE_SStructBiCGSTABDestroy(solver);

		// save the result into the domain objects
		dc.saveResult();
		double ausum2 = dc.integrateAU();
		double fsum2  = dc.integrateF();
		double bflux  = dc.integrateBoundaryFlux();
		if (my_global_rank == 0) {
			std::cout << u8"Σf-Au: " << fsum2 - ausum2 << endl;
			std::cout << u8"Σf: " << fsum2 << endl;
			std::cout << u8"ΣAu: " << ausum2 << endl;
			if (f_neumann) {
				std::cout << u8"∮ du/dn: " << bflux << endl;
				std::cout << u8"∮ du/dn - Σf: " << bflux - fsum2 << endl;
				std::cout << u8"∮ du/dn - ΣAu: " << bflux - ausum2 << endl;
			}
		}

		// Calcuate error
		double exact_norm;
		double diff_norm;

		if (f_neumann) {
			double uavg = dc.integrateU()/dc.area();
			double eavg = dc.integrateExact()/dc.area();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			exact_norm = dc.exactNorm(eavg);
			diff_norm  = dc.diffNorm(uavg, eavg);
		} else {
			exact_norm = dc.exactNorm();
			diff_norm  = dc.diffNorm();
		}

		MPI_Barrier(MPI_COMM_WORLD);
		steady_clock::time_point total_stop = steady_clock::now();
		duration<double>         total_time = total_stop - domain_start;

		double residual = dc.residual();
		double fnorm    = dc.fNorm();
		double ausum    = dc.integrateAU();
		double fsum     = dc.integrateF();
		if (my_global_rank == 0) {
			std::cout << "Total run time: " << total_time.count() << endl;
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << diff_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum-fsum << endl;
			// if (f_neumann) {
			//	std::cout << u8"∮ du/dn - ΣAu: " << dc.sumBoundaryFlux() - ausum << endl;
			//}
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
            //TODO
		}
        if (f_flux){
            //TODO
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
		cout << endl;
		cout << "Average: " << times.sum() / times.size() << endl;
	}

	return 0;
}
