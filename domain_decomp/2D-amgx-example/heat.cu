#include "DomainSignatureCollection.h"
#include "DomainCollection.h"
#include "Timer.h"
#include "args.h"
#include "amgx_c.h"
#include <chrono>
#include <cmath>
#include <iostream>
//#include <mpi.h>
#include <string>
#include <fstream>
#include <unistd.h>
#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#endif

// =========== //
// main driver //
// =========== //

using namespace std;
void print_callback(const char *msg, int length)
{
    cout << msg;
}
int main(int argc, char *argv[])
{
	using namespace std::chrono;


	MPI_Init(NULL, NULL);
	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// parse input
	args::ArgumentParser parser("");
	args::HelpFlag       help(parser, "help", "Display this help menu", {'h', "help"});

	args::ValueFlag<int> f_n(parser, "n", "number of cells in the x direction, in each domain",
	                         {'n'});
	args::ValueFlag<string> f_mesh(parser, "file_name", "read in a mesh", {"mesh"});
	args::ValueFlag<int>    f_square(parser, "num_domains",
	                              "create a num_domains x num_domains square of grids", {"square"});
	args::ValueFlag<int> f_amr(parser, "num_domains", "create a num_domains x num_domains square "
	                                                  "of grids, and a num_domains*2 x "
	                                                  "num_domains*2 refined square next to it",
	                           {"amr"});
	args::Flag f_outclaw(parser, "outclaw", "output amrclaw ascii file", {"outclaw"});
#ifdef HAVE_VTK
	args::Flag f_outvtk(parser, "", "output to vtk format", {"outvtk"});
#endif
	args::ValueFlag<int> f_l(parser, "n", "run the program n times and print out the average",
	                         {'l'});
	args::ValueFlag<string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                            {'m'});
	args::ValueFlag<string> f_s(parser, "solution filename", "the file to write the solution to",
	                            {'s'});
	args::ValueFlag<string> f_resid(parser, "residual filename",
	                                "the file to write the residual to", {"residual"});
	args::ValueFlag<string> f_error(parser, "error filename", "the file to write the error to",
	                                {"error"});
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
	args::ValueFlag<double> f_omega(
	parser, "tolerance", "set the tolerance of the iterative solver (default is 1e-10)", {"omega"});
	args::ValueFlag<int> f_d(
	parser, "row", "pin gamma value to zero (by modifying that row of the schur compliment matrix)",
	{'z'});
	args::ValueFlag<int> f_div(parser, "divide", "use iterative method", {"divide"});
	args::Flag           f_wrapper(parser, "wrapper", "use a function wrapper", {"wrap"});
	args::Flag           f_blockcrs(parser, "wrapper", "use a function wrapper", {"blockcrs"});
	args::Flag           f_crs(parser, "wrapper", "use a function wrapper", {"crs"});
	args::Flag           f_gauss(parser, "gauss", "solve gaussian function", {"gauss"});
	args::Flag           f_zero(parser, "gauss", "solve gaussian function", {"zero"});
	args::Flag           f_prec(parser, "prec", "use block diagonal preconditioner", {"prec"});
	args::Flag           f_precblockj(parser, "prec", "use block diagonal jacobi preconditioner",
	                        {"precblockj"});
	args::Flag f_precj(parser, "prec", "use block diagonal jacobi preconditioner", {"precj"});
	args::Flag f_precmuelu(parser, "prec", "use AMG preconditioner", {"muelu"});
	args::Flag f_precddmg(parser, "prec", "use AMG preconditioner", {"ddmg"});
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {"neumann"});
	args::Flag f_cg(parser, "gmres", "use CG for iterative solver", {"cg"});
	args::Flag f_gmres(parser, "gmres", "use GMRES for iterative solver", {"gmres"});
	args::Flag f_lsqr(parser, "gmres", "use GMRES for iterative solver", {"lsqr"});
	args::Flag f_rgmres(parser, "rgmres", "use GCRO-DR (Recycling GMRES) for iterative solver",
	                    {"rgmres"});
	args::Flag f_bicg(parser, "gmres", "use BiCGStab for iterative solver", {"bicg"});
	args::Flag f_zerou(parser, "zerou", "modify matrix so that the sum of the solution is zero",
	                   {"nozerou"});
	args::Flag f_nozerof(parser, "zerou", "don't  make the rhs match the boundary conditions",
	                     {"nozerof"});
	args::Flag f_pingamma(parser, "pingamma", "pin the first gamma to zero", {"pingamma"});
	args::Flag f_lu(parser, "lu", "use KLU solver", {"klu"});
	args::Flag f_mumps(parser, "lu", "use MUMPS solver", {"mumps"});
	args::Flag f_basker(parser, "lu", "use Basker solver", {"basker"});
	args::Flag f_superlu(parser, "lu", "use SUPERLU solver", {"superlu"});
	args::Flag f_ilu(parser, "ilu", "use incomplete LU preconditioner", {"ilu"});
	args::Flag f_riluk(parser, "ilu", "use RILUK preconditioner", {"riluk"});
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

	bool direct_solve = (f_lu || f_superlu || f_mumps || f_basker);
	bool use_crs = (f_crs || direct_solve || f_ilu || f_riluk || f_precj || f_precmuelu || f_prec);

	DomainSignatureCollection dsc;
	if (f_mesh) {
		string d = args::get(f_mesh);
		dsc      = DomainSignatureCollection(d, my_global_rank);
	} else if (f_amr) {
		int d = args::get(f_amr);
		dsc   = DomainSignatureCollection(d, d, my_global_rank, true);
	} else {
		int d = args::get(f_square);
		dsc   = DomainSignatureCollection(d, d, my_global_rank);
	}
	if (f_div) {
		for (int i = 0; i < args::get(f_div); i++) {
			dsc.divide();
		}
	}
	if (f_neumann) {
		dsc.setNeumann();
		if (!f_pingamma && !f_zerou) {
		dsc.setZeroPatch();
		}
	}
	// Set the number of discretization points in the x and y direction.
	int nx          = args::get(f_n);
	int ny          = args::get(f_n);
	int total_cells = dsc.num_global_domains * nx * ny;
	cerr << "Total cells: " << total_cells << endl;

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

	if (f_zero) {
		ffun  = [](double x, double y) { return 0; };
		gfun  = [](double x, double y) { return 0; };
		nfunx = [](double x, double y) { return 0; };
		nfuny = [](double x, double y) { return 0; };
	} else if (f_gauss) {
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

    //library handles
    AMGX_Mode mode;
    AMGX_config_handle cfg;
    AMGX_resources_handle rsrc;
    AMGX_matrix_handle gA;
    AMGX_vector_handle gb, gx;
    AMGX_solver_handle solver;
    mode = AMGX_mode_dDDI;
    //status handling
    AMGX_SOLVE_STATUS status;
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());
    /* create resources, matrix, vector and solver */
    AMGX_config_create_from_file(&cfg, "amgx.json");
    AMGX_resources_create_simple(&rsrc, cfg);
    AMGX_matrix_create(&gA, rsrc, mode);
    AMGX_vector_create(&gx, rsrc, mode);
    AMGX_vector_create(&gb, rsrc, mode);
    AMGX_solver_create(&solver, rsrc, mode, cfg);

	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		DomainCollection dc(dsc, nx);
		if (f_neumann) {
			if (f_neumann && f_zerou) {
				dc.setZeroU();
			}
		}

		if (f_neumann) {
			dc.initNeumann(ffun, gfun, nfunx, nfuny, f_amr);
		} else {
			dc.initDirichlet(ffun, gfun);
			dc.amr = f_amr;
		}

		timer.stop("Domain Initialization");

		// Create a map that will be used in the iterative solver

        int num_rows = nx*dsc.iface_map_vec.size();
		// Create the gamma and diff vectors
		double*                   gamma = new double[num_rows];
            for(int i=0;i<num_rows;i++){
                gamma[i]=0;
            }
		double*                   r     = new double[num_rows];
		double*                   x     = new double[num_rows];
		double*                   d     = new double[num_rows];
		double*                   diff  = new double[num_rows];
		double*                   b     = new double[num_rows];
		AmgxCrs                   A;

		if (f_neumann && !f_nozerof) {
			double fdiff = (dc.integrateBoundaryFlux() - dc.integrateF()) / dc.area();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
            dc.zeroF(fdiff);
		}
		timer.start("Complete Solve");
		if (dsc.num_global_domains != 1) {
			// do iterative solve

			// Get the b vector
			dc.solveWithInterface(gamma, b);

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

            timer.start("Matrix Formation");

            A= dc.formCRSMatrix();

            timer.stop("Matrix Formation");
            timer.start("AMGX Setup");

            
            AMGX_SAFE_CALL(AMGX_vector_upload(gb,num_rows,1,(void*)b));
            AMGX_SAFE_CALL(AMGX_vector_set_zero(gx,num_rows,1));
            AMGX_SAFE_CALL(AMGX_matrix_upload_all(gA,num_rows,A.nnz,1,1,&A.row_ptrs[0],&A.cols[0],(void*)&A.data[0],nullptr));
            AMGX_SAFE_CALL(AMGX_solver_setup(solver, gA));

            timer.stop("AMGX Setup");

			timer.stop("Linear System Setup");

		    timer.start("Gamma Solve");
            AMGX_SAFE_CALL(AMGX_solver_solve_with_0_initial_guess(solver,gb,gx));
            AMGX_SAFE_CALL(AMGX_vector_download(gx,(void*)gamma));

		    timer.stop("Gamma Solve");
		}

		// Do one last solve
		timer.start("Patch Solve");

		dc.solveWithInterface(gamma, diff);

		timer.stop("Patch Solve");

		dc.residual();
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
        /*
           TODO iterative
		if (f_iter && !direct_solve) {
			timer.start("Iterative Refinement Step");
			dc.residual();
			dc.swapResidSol();

			if (dsc.num_global_domains != 1) {
				x->putScalar(0);
				dc.solveWithInterface(*x, *r);
				// op->apply(*gamma, *r);
				// r->update(1.0, *b, -1.0);

				solver->reset(Belos::ResetType::Problem);
				if (f_wrapper) {
					((FuncWrap *) op.getRawPtr())->setB(r);
				}
				problem->setProblem(x, r);
				solver->setProblem(problem);
				solver->solve();
			}
			dc.solveWithInterface(*x, *d);
			dc.sumResidIntoSol();
			timer.stop("Iterative Refinement Step");
		}
        */

		///////////////////
		// solve end
		///////////////////
		timer.stop("Complete Solve");

		// Calcuate error
        double exact_norm;
        double diff_norm;

		if (f_neumann) {
			double uavg = dc.integrateU() / dc.area();
			double eavg = dc.integrateExact() / dc.area();

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

		double residual = dc.residual();
		double fnorm    = dc.fNorm();
		double ausum    = dc.integrateAU();
		double fsum     = dc.integrateF();
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << diff_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			// if (f_neumann) {
			//	std::cout << u8"∮ du/dn - ΣAu: " << dc.sumBoundaryFlux() - ausum << endl;
			//}
			cout.unsetf(std::ios_base::floatfield);
		}
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
#ifdef HAVE_VTK
		if (f_outvtk) {
			dc.outputVTK();
		}
#endif
		cout.unsetf(std::ios_base::floatfield);
	}

	if (my_global_rank == 0) {
		cout << timer;
	}
	return 0;
}
