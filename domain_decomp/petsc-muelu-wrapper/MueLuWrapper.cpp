#include "MueLuWrapper.h"
#include "MyTypeDefs.h"
#include "PW.h"
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_TpetraOperator_def.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <map>
#include <set>
#include <vector>
using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
struct Vars {
	RCP<matrix_type>              A;
	RCP<vec_type>                 x;
	RCP<vec_type>                 b;
	RCP<map_type>                 map;
	RCP<const Teuchos::Comm<int>> comm;
	RCP<Belos::LinearProblem<scalar_type, vec_type, Tpetra::Operator<scalar_type>>> problem;
	RCP<Belos::SolverManager<scalar_type, vec_type, Tpetra::Operator<scalar_type>>> solver;
};
MueLuWrapper::MueLuWrapper(Mat A, double tol, string config)
{
	vars       = new Vars;
	vars->comm = rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	PW<Mat> localA;
	MatMPIAIJGetLocalMat(A, MAT_INITIAL_MATRIX, &localA);

	int m, n;
	MatGetLocalSize(A, &m, &n);
	num_rows = n;

	const int *ranges;
	MatGetOwnershipRanges(A, &ranges);
	vector<int> procs;
	for (int p = 0; p < size; p++) {
		for (int i = ranges[n]; i < ranges[n + 1]; i++) {
			procs.push_back(p);
		}
	}

	// create map
	vars->map = rcp(new map_type(-1, num_rows, 0, vars->comm));

	// create matrix
	vars->A = rcp(new matrix_type(vars->map, 10));

	for (int i = 0; i < num_rows; i++) {
		int           row = i + ranges[rank];
		int           ncols;
		const int *   c;
		const double *vals;
		MatGetRow(localA, i, &ncols, &c, &vals);
		vars->A->insertGlobalValues(row, ncols, vals, c);
	}

	vars->A->fillComplete();

	// creat vectors
	vars->x = rcp(new vec_type(vars->map, 1));
	vars->b = rcp(new vec_type(vars->map, 1));

	// create muelu preconditioner
	vars->problem
	= rcp(new Belos::LinearProblem<scalar_type, vec_type, Tpetra::Operator<scalar_type>>(
	vars->A, vars->x, vars->b));
	vars->problem->setLeftPrec(MueLu::CreateTpetraPreconditioner(vars->A, config));

	Teuchos::ParameterList belosList;
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
	vars->solver
	= rcp(new Belos::BlockGmresSolMgr<scalar_type, vec_type, Tpetra::Operator<scalar_type>>(
	vars->problem, rcp(&belosList, false)));
}
MueLuWrapper::~MueLuWrapper() { delete vars; }
void MueLuWrapper::solve(Vec x, Vec b)
{
	double *x_ptr, *b_ptr;
	VecGetArray(x, &x_ptr);
	VecGetArray(b, &b_ptr);
	auto tx_view = vars->x->get1dViewNonConst();
	auto tb_view = vars->b->get1dViewNonConst();
	for (int i = 0; i < num_rows; i++) {
		tx_view[i] = x_ptr[i];
		tb_view[i] = b_ptr[i];
	}
	vars->problem->setProblem();
	vars->solver->solve();
	for (int i = 0; i < num_rows; i++) {
		x_ptr[i] = tx_view[i];
	}
	VecRestoreArray(x, &x_ptr);
	VecRestoreArray(b, &b_ptr);
}
