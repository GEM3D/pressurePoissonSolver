#include "DomainCollection.h"
#include "Factory.h"
#include "SchurHelper.h"
//#include "FunctionWrapper.h"
#include "FivePtPatchOperator.h"
#include "Init.h"
#include "MyTypeDefs.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "PatchSolvers/FishpackPatchSolver.h"
#include "QuadInterpolator.h"
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
#include <Ifpack2_BlockRelaxation_decl.hpp>
#include <Ifpack2_BlockRelaxation_def.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_ILUT_decl.hpp>
#include <Ifpack2_ILUT_def.hpp>
#include <Ifpack2_Relaxation_decl.hpp>
#include <Ifpack2_Relaxation_def.hpp>
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
#include <cmath>
#include <iostream>
#include <memory>
#include <petscvec.h>
#include <string>
#include <unistd.h>
#include <petscksp.h>
#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#endif

using Teuchos::RCP;
using Teuchos::rcp;

// =========== //
// main driver //
// =========== //

using namespace std;

int main(int argc, char *argv[])
{
    PetscInitialize(&argc,&argv,nullptr,nullptr);
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
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {"neumann"});

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

	// preconditioners
	args::Flag f_precj(parser, "prec", "use jacobi preconditioner", {"precj"});
	args::Flag f_prec(parser, "prec", "use block diagonal jacobi preconditioner", {"prec"});
	args::Flag f_precmuelu(parser, "prec", "use MueLu AMG preconditioner", {"muelu"});
	// iterative options
	args::Flag f_cg(parser, "gmres", "use CG for iterative solver", {"cg"});
	args::Flag f_gmres(parser, "gmres", "use GMRES for iterative solver", {"gmres"});
	args::Flag f_lsqr(parser, "gmres", "use least squares for iterative solver", {"lsqr"});
	args::Flag f_rgmres(parser, "rgmres", "use GCRO-DR (Recycling GMRES) for iterative solver",
	                    {"rgmres"});
	args::Flag f_bicg(parser, "gmres", "use BiCGStab for iterative solver", {"bicg"});

	// direct solvers
	args::Flag f_lu(parser, "lu", "use KLU solver", {"klu"});
	args::Flag f_mumps(parser, "lu", "use MUMPS solver", {"mumps"});
	args::Flag f_basker(parser, "lu", "use Basker solver", {"basker"});
	args::Flag f_superlu(parser, "lu", "use SUPERLU solver", {"superlu"});
	args::Flag f_ilu(parser, "ilu", "use incomplete LU preconditioner", {"ilu"});
	args::Flag f_riluk(parser, "ilu", "use RILUK preconditioner", {"riluk"});
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
	// Set the number of discretization points in the x and y direction.
	int nx = args::get(f_n);
	int ny = args::get(f_n);

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

	bool direct_solve = (f_lu || f_superlu || f_mumps || f_basker);

	///////////////
	// Create Mesh
	///////////////
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
	if (f_neumann) {
		dc.setNeumann();
	}
	dc.n = nx;
	for (auto &p : dc.domains) {
		p.second.n = nx;
	}
	int total_cells = dc.getGlobalNumCells();
	cout << "Total cells: " << total_cells << endl;

	if (dc.num_global_domains < num_procs) {
		std::cerr << "number of domains must be greater than or equal to the number of processes\n";
		return 1;
	}
	// partition domains if running in parallel
	if (num_procs > 1) {
		dc.zoltanBalance();
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
		ffun
		= [](double x, double y) { return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x); };
		gfun  = [](double x, double y) { return sinl(M_PI * y) * cosl(2 * M_PI * x); };
		nfunx = [](double x, double y) { return -2 * M_PI * sinl(M_PI * y) * sinl(2 * M_PI * x); };
		nfuny = [](double x, double y) { return M_PI * cosl(M_PI * y) * cosl(2 * M_PI * x); };
	}

	// set the patch solver
	RCP<PatchSolver> p_solver;
	if (f_fish) {
		p_solver = rcp(new FishpackPatchSolver());
	} else {
		p_solver = rcp(new FftwPatchSolver(dc));
	}

	// patch operator
	RCP<PatchOperator> p_operator = rcp(new FivePtPatchOperator());

	// interface interpolator
	RCP<Interpolator> p_interp = rcp(new QuadInterpolator());

	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		SchurHelper sch(dc, comm, p_solver, p_operator, p_interp);

		IS              domain_is = dc.getDomainIS();
	    Vec u         = dc.getNewDomainVec();
	    Vec exact     = dc.getNewDomainVec();
	    Vec f         = dc.getNewDomainVec();

		if (f_neumann) {
			Init::initNeumann(dc, nx, f, exact, ffun, gfun, nfunx, nfuny);
		} else {
			Init::initDirichlet(dc, nx, f, exact, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		// Create the gamma and diff vectors
		shared_ptr<Vec>                    gamma = dc.getNewSchurVec();
		shared_ptr<Vec>                    diff  = dc.getNewSchurVec();
		shared_ptr<Vec>                    b     = dc.getNewSchurVec();
		shared_ptr<Mat>                   A;

		// Create linear problem for the Belos solver
        shared_ptr<KSP>solver (new KSP,KSPDestroy);
        //KSPCreate(MPI_COMM_WORLD,solver.get());

		if (f_neumann && !f_nozerof) {
            /*
			double fdiff = (dc.integrate(*f)) / dc.area();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;

			RCP<vector_type> diff = rcp(new vector_type(domain_map, 1));
			diff->putScalar(fdiff);
			f->update(1.0, *diff, 1.0);
            */
		}

#ifdef ENABLE_AMGX
		Teuchos::RCP<AmgxWrapper> amgxsolver;
#endif
#ifdef ENABLE_HYPRE
		Teuchos::RCP<HypreWrapper> hypresolver;
#endif

		if (dc.num_global_domains != 1) {
			// do iterative solve

			// Get the b vector
			sch.solveWithInterface(f, u, *gamma, *b);
            VecScale(*b,-1.0);

			if (f_rhs) {
		//		Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(args::get(f_rhs), b, "",
		//		                                                          "");
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

				A = sch.formCRSMatrix();

				timer.stop("Matrix Formation");

				if (f_m) {
					//Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(args::get(f_m), A,
					 //                                                          "", "");
				}
			}

			if (f_amgx) {
#ifdef ENABLE_AMGX
				timer.start("AMGX Setup");
		//		amgxsolver = rcp(new AmgxWrapper(A, dc, nx));
				timer.stop("AMGX Setup");
#endif
			} else if (f_hypre) {
#ifdef ENABLE_HYPRE
				timer.start("Hypre Setup");
		//		hypresolver = rcp(new HypreWrapper(A, dc, nx, tol, true));
				timer.stop("Hypre Setup");
#endif
			} else {
                //preconditoners
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

		if (dc.num_global_domains != 1) {
			timer.start("Gamma Solve");
			if (f_amgx) {
// solve
#ifdef ENABLE_AMGX
				//amgxsolver->solve(gamma, b);
#endif
			} else if (f_hypre) {
#ifdef ENABLE_HYPRE
				//hypresolver->solve(gamma, b);
#endif
			} else {

                KSPSetTolerances(*solver,tol,PETSC_DEFAULT,PETSC_DEFAULT,5000);
                KSPSetOperators(*solver,*A,*A);
                KSPSolve(*solver,*b,*gamma);
			}
			timer.stop("Gamma Solve");
		}

		// Do one last solve
		timer.start("Patch Solve");

		sch.solveWithInterface(f, u, *gamma, *diff);

		timer.stop("Patch Solve");

		///////////////////
		// solve end
		///////////////////
		timer.stop("Complete Solve");

		// residual
        Vec resid = dc.getNewDomainVec();
        Vec au = dc.getNewDomainVec();
		sch.applyWithInterface(u, *gamma, au);
        VecAXPBYPCZ(resid,-1.0,1.0,0.0,au,f);
		double residual;
        VecNorm(resid,NORM_2,&residual);
		double fnorm;
        VecNorm(f,NORM_2,&fnorm);

		// error
        Vec error = dc.getNewDomainVec();
        VecAXPBYPCZ(error,-1.0,1.0,0.0,exact,u);
		if (f_neumann) {
		/*	double uavg = dc.integrate(*u) / dc.area();
			double eavg = dc.integrate(*exact) / dc.area();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			vector_type ones(domain_map, 1);
			ones.putScalar(1);
			error->update(eavg - uavg, ones, 1.0);*/
		}
		double error_norm;
        VecNorm(error,NORM_2,&error_norm);
		double exact_norm;
        VecNorm(exact,NORM_2,&exact_norm);

		double ausum = dc.integrate(au);
		double fsum  = dc.integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
		}

		// output
		MMWriter mmwriter(dc, f_amr);
		if (f_s) {
			//mmwriter.write(*u, args::get(f_s));
		}
		if (f_resid) {
			//mmwriter.write(*resid, args::get(f_resid));
		}
		if (f_error) {
			//mmwriter.write(*error, args::get(f_error));
		}
		if (f_g) {
			//Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(args::get(f_g), gamma, "",
			 //                                                         "");
		}
		if (f_outclaw) {
			ClawWriter writer(dc);
			//writer.write(*u, *resid);
		}
#ifdef HAVE_VTK
		if (f_outvtk) {
			VtkWriter writer(dc, args::get(f_outvtk));
			//writer.add(*u, "Solution");
			//writer.add(*error, "Error");
			//writer.add(*resid, "Residual");
			writer.write();
		}
#endif
		cout.unsetf(std::ios_base::floatfield);
	}

	if (my_global_rank == 0) {
		cout << timer;
	}
	return 0;
}
