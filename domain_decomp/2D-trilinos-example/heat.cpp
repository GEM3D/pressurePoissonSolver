#include "BlockJacobiRelaxer.h"
#include "DomainSignatureCollection.h"
#include "Factory.h"
#include "FunctionWrapper.h"
#include "MyTypeDefs.h"
#include "OpShift.h"
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
#include <chrono>
#include <cmath>
#include <iostream>
#include <iostream>
//#include <mpi.h>
#include <string>
#include <unistd.h>
#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#endif

#include "DDMultiGrid.h"

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
	args::Flag           f_outclaw(parser, "outclaw", "output amrclaw ascii file", {"outclaw"});
	args::Flag           f_outvtk(parser, "", "output to vtk format", {"outvtk"});
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
		dsc      = DomainSignatureCollection(d, comm->getRank());
	} else if (f_amr) {
		int d = args::get(f_amr);
		dsc   = DomainSignatureCollection(d, d, comm->getRank(), true);
	} else {
		int d = args::get(f_square);
		dsc   = DomainSignatureCollection(d, d, comm->getRank());
	}
	if (f_div) {
		for (int i = 0; i < args::get(f_div); i++) {
			dsc.divide();
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

	scalar_type tol = 1e-10;
	if (f_t) {
		tol = args::get(f_t);
	}

	if (f_omega) {
		BlockJacobiRelaxer::omega = args::get(f_omega);
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

	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		DomainCollection dc(dsc, nx, comm);
		if (f_neumann) {
			if (f_neumann && f_zerou) {
				dc.setZeroU();
			} else if (!f_pingamma && dc.domains.count(0)) {
				dc.domains[0]->zero_patch = true;
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
		RCP<map_type>       matrix_map       = dc.matrix_map;
		RCP<const map_type> matrix_map_const = dc.matrix_map;

		// Create the gamma and diff vectors
		RCP<vector_type>                   gamma = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   r     = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   x     = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   d     = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   diff  = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   b     = rcp(new vector_type(matrix_map, 1));
		RCP<RBMatrix>                      RBA;
		RCP<matrix_type>                   A;
		RCP<single_vector_type>            s;
		RCP<Tpetra::Operator<scalar_type>> op;
		RCP<const Tpetra::RowMatrix<>>     rm;
		RCP<Amesos2::Solver<matrix_type, vector_type>> dsolver;

		// Create linear problem for the Belos solver
		RCP<Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>> problem;
		RCP<Belos::SolverManager<scalar_type, vector_type, Tpetra::Operator<scalar_type>>> solver;
		Teuchos::ParameterList belosList;

		typedef Ifpack2::Preconditioner<scalar_type> Preconditioner;
		RCP<Preconditioner>                          prec;
		if (f_neumann && !f_nozerof) {
			double fdiff = (dc.integrateBoundaryFlux() - dc.integrateF()) / dc.area();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;
			for (auto &p : dc.domains) {
				Domain &d = *p.second;
				d.f += fdiff;
			}
		}
		steady_clock::time_point tsolve_start;
		if (dsc.num_global_domains != 1) {
			// do iterative solve

			// Get the b vector
			dc.solveWithInterface(*gamma, *b);

			if (save_rhs_file != "") {
				Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_rhs_file, b, "", "");
			}

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			if (use_crs) {
				timer.start("Matrix Formation");

				dc.formCRSMatrix(matrix_map, A, &s);

				timer.stop("Matrix Formation");

				if (f_neumann && f_zerou) {
					RCP<OpShift> os = rcp(new OpShift(A, s));
					op              = os;
				} else {
					op = A;
				}
				rm = A;

				if (save_matrix_file != "")
					Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(save_matrix_file, A,
					                                                           "", "");
			} else if (f_wrapper) {
				// Create a function wrapper
				op = rcp(new FuncWrap(b, &dc));
			} else {
				// Form the matrix
				timer.start("Matrix Formation");

				dc.formRBMatrix(matrix_map, RBA, &s);

				timer.stop("Matrix Formation");

				if (save_matrix_file != "") {
					comm->barrier();
					steady_clock::time_point write_start = steady_clock::now();

					ofstream out_file(save_matrix_file);
					out_file << *RBA;
					out_file.close();

					// Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(
					// save_matrix_file + ".s", RBA->shift_vec, "", "");

					comm->barrier();
					duration<double> write_time = steady_clock::now() - write_start;
					if (my_global_rank == 0)
						cout << "Time to write matrix to file: " << write_time.count() << endl;
				}

				if (f_neumann && f_zerou) {
					RCP<OpShift> os = rcp(new OpShift(RBA, s));
					op              = os;
				} else {
					op = RBA;
				}
			}

			problem
			= rcp(new Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>(
			op, gamma, b));

			if (f_precmuelu) {
				timer.start("MueLu Preconditioner Formation");

				Teuchos::RCP<vector_type> xy = dc.getInterfaceCoords();
				Teuchos::RCP<op_type>     P  = Factory::getAmgPreconditioner(A, xy);

				problem->setLeftPrec(P);

				timer.stop("MueLu Preconditioner Formation");

			} else if (f_riluk) {
				timer.start("RILUK Preconditioner Formation");

				Teuchos::RCP<Ifpack2::RILUK<Tpetra::RowMatrix<>>> P
				= rcp(new Ifpack2::RILUK<Tpetra::RowMatrix<>>(rm));
				P->compute();

				problem->setLeftPrec(P);

				timer.stop("RILUK Preconditioner Formation");
			} else if (f_ilu) {
				timer.start("ILUT Preconditioner Formation");

				Teuchos::RCP<Ifpack2::ILUT<Tpetra::RowMatrix<>>> P
				= rcp(new Ifpack2::ILUT<Tpetra::RowMatrix<>>(rm));
				Teuchos::ParameterList params;
				params.set("fact: ilut level-of-fill", 3);
				params.set("fact: drop tolerance", 0.0);
				params.set("fact: absolute threshold", 0.1);
				// P->setParameters(params);

				P->compute();
				problem->setLeftPrec(P);

				timer.stop("ILUT Preconditioner Formation");
			} else if (f_precj) {
				timer.start("Jacobi Preconditioner Formation");

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

				timer.stop("Jacobi Preconditioner Formation");
			} else if (f_precblockj) {
				timer.start("Block Jacobi Preconditioner Formation");

				RCP<BlockJacobiRelaxer> P = rcp(new BlockJacobiRelaxer(RBA, s));
				problem->setLeftPrec(P);

				timer.stop("Block Jacobi Preconditioner Formation");
			} else if (f_precddmg) {
				timer.start("DDMG Preconditioner Formation");

				RCP<DDMultiGrid> P = rcp(new DDMultiGrid(&dc, RBA));

				problem->setRightPrec(P);

				timer.stop("DDMG Preconditioner Formation");
			} else if (f_prec) {
				timer.start("Block Diagonal Preconditioner Formation");

				// Create the relaxation.  You could also do this using
				// Ifpack2::Factory (the preconditioner factory) if you like.
				RCP<precond_type> prec = rcp(new Ifpack2::BlockRelaxation<Tpetra::RowMatrix<>>(A));
				// Make the list of relaxation parameters.
				Teuchos::ParameterList params;
				// Do symmetric SOR / Gauss-Seidel.
				params.set("relaxation: type", "Jacobi");
				// Two sweeps (of symmetric SOR / Gauss-Seidel) per apply() call.
				params.set("relaxation: sweeps", 1);
				params.set("relaxation: container", "Dense");
				// ... Set any other parameters you want to set ...
				params.set("partitioner: local parts", dsc.matrix_j_high - dsc.matrix_j_low);

				// Set parameters.
				prec->setParameters(params);
				// Prepare the relaxation instance for use.
				prec->initialize();
				prec->compute();

				timer.stop("Block Diagonal Preconditioner Formation");
			}

			if (direct_solve) {
				timer.start("LU Factorization");
				string name;
				if (f_lu) {
					name = "KLU2";
				}
				if (f_mumps) {
					name = "mumps";
				}
				if (f_superlu) {
					name = "superlu_dist";
				}
				if (f_basker) {
					name = "Basker";
				}

				dsolver = Amesos2::create<matrix_type, vector_type>(name, A);

				dsolver->symbolicFactorization().numericFactorization();
				timer.stop("LU Factorization");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");

			///////////////////
			// solve start
			///////////////////
			timer.start("Complete Solve");

			timer.start("Gamma Solve");
			if (f_read_gamma) {
				gamma = Tpetra::MatrixMarket::Reader<matrix_type>::readDenseFile(
				args::get(f_read_gamma), comm, matrix_map_const);
			} else if (f_lu || f_superlu || f_mumps || f_basker) {
				dsolver->solve(&*gamma, &*b);
				if (f_iter) {
					x->putScalar(0);
					A->apply(*gamma, *r);
					r->update(1.0, *b, -1.0);

					dsolver->solve(&*x, &*r);

					gamma->update(1.0, *x, 1.0);
				}

			} else {
				problem->setProblem();

				// Set the parameters
				belosList.set("Block Size", 1);
				belosList.set("Maximum Iterations", 5000);
				belosList.set("Convergence Tolerance", tol);
				belosList.set("Output Frequency", 1);
				int verbosity = Belos::Errors + /* Belos::StatusTestDetails +*/ Belos::Warnings
				                + Belos::TimingDetails + Belos::Debug + Belos::IterationDetails;
				belosList.set("Verbosity", verbosity);
				// belosList.set("Orthogonalization", "ICGS");
				//
				belosList.set("Rel RHS Err", 0.0);
				belosList.set("Rel Mat Err", 0.0);

				// Create solver and solve
				if (f_rgmres) {
					solver = rcp(new Belos::GCRODRSolMgr<scalar_type, vector_type,
					                                     Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				} else if (f_cg) {
					solver = rcp(new Belos::BlockCGSolMgr<scalar_type, vector_type,
					                                      Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				} else if (f_bicg) {
					solver = rcp(new Belos::BiCGStabSolMgr<scalar_type, vector_type,
					                                       Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				} else {
					solver = rcp(new Belos::BlockGmresSolMgr<scalar_type, vector_type,
					                                         Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				}
				solver->solve();
			}
		}
		timer.stop("Gamma Solve");

		// Do one last solve
		timer.start("Patch Solve");

		dc.solveWithInterface(*gamma, *diff);

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

		///////////////////
		// solve end
		///////////////////
		timer.stop("Complete Solve");

		// Calcuate error
		RCP<map_type>    err_map = rcp(new map_type(-1, 1, 0, comm));
		Tpetra::Vector<> exact_norm(err_map);
		Tpetra::Vector<> diff_norm(err_map);

		if (f_neumann) {
			double uavg = dc.integrateU() / dc.area();
			double eavg = dc.integrateExact() / dc.area();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
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

		double residual = dc.residual();
		double fnorm    = dc.fNorm();
		double ausum    = dc.integrateAU();
		double fsum     = dc.integrateF();
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << global_diff_norm / global_exact_norm << endl;
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
		if (save_gamma_file != "") {
			Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_gamma_file, gamma, "",
			                                                          "");
		}
		if (f_flux) {
			dc.getFluxDiff(*gamma);
			Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(args::get(f_flux), gamma, "",
			                                                          "");
		}
		if (f_outclaw) {
			dc.outputClaw();
		}
		if (f_outvtk) {
			dc.outputVTK();
		}
		cout.unsetf(std::ios_base::floatfield);
	}

	if (my_global_rank == 0) {
		cout << timer;
	}
	return 0;
}
