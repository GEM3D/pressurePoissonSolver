#include "args.h"
#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_CombineMode.h>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <string>
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>


using Teuchos::RCP;
using Teuchos::rcp;
typedef Epetra_CrsMatrix   matrix_type;
typedef Epetra_MultiVector vector_type;
typedef Amesos2::Solver<matrix_type, vector_type> solver_type;
typedef Epetra_Map map_type;
class Domain
{
	public:
	RCP<matrix_type> A;
	RCP<vector_type> f;
	RCP<vector_type> u;
	RCP<solver_type> solver;
	int              nx;
	int              ny;
	double           h_x;
	double           h_y;
	double *         boundary_north = nullptr;
	double *         boundary_south = nullptr;
	double *         boundary_east  = nullptr;
	double *         boundary_west  = nullptr;
	bool             has_north      = false;
	bool             has_south      = false;
	bool             has_east       = false;
	bool             has_west       = false;
	RCP<map_type>    domain_map;

	Domain(RCP<matrix_type> A, RCP<vector_type> f, int nx, int ny, double h_x, double h_y)
	{
		solver = Amesos2::create<matrix_type, vector_type>("KLU2", A);
		solver->symbolicFactorization().numericFactorization();
		this->A   = A;
		this->f   = f;
		this->nx  = nx;
		this->ny  = ny;
		this->h_x = h_x;
		this->h_y = h_y;
		u         = rcp(new vector_type(f->Map(), 1, false));
	}

	void solve()
	{
		RCP<vector_type> f_copy = rcp(new vector_type(f->Map(), 1, false));
		for (int x_i = 0; x_i < nx; x_i++) {
			for (int y_i = 0; y_i < ny; y_i++) {
				(*f_copy)[0][y_i * nx + x_i] = (*f)[0][y_i * nx + x_i];
				if (y_i == ny - 1 && has_north) {
					(*f_copy)[0][y_i * nx + x_i] += -2.0 / (h_y * h_y) * boundary_north[x_i];
				}
				if (x_i == nx-1 && has_east) {
					(*f_copy)[0][y_i * nx + x_i] += -2.0 / (h_x * h_x) * boundary_east[y_i];
				}
				if (y_i == 0 && has_south) {
					(*f_copy)[0][y_i * nx + x_i] += -2.0 / (h_y * h_y) * boundary_south[x_i];
				}
				if (x_i == 0 && has_west) {
					(*f_copy)[0][y_i * nx + x_i] += -2.0 / (h_x * h_x) * boundary_west[y_i];
				}
			}
		}
		solver->setX(u);
		solver->setB(f_copy);
		solver->solve();
	}

	void solveWithInterface(const vector_type &gamma, vector_type &diff)
	{
		int my_global_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);
		if (my_global_rank == 0) {
			std::cerr << "i have been called!\n";
		}
		Epetra_Import importer(*domain_map, gamma.Map());
		vector_type   local_vector(*domain_map, 1);
		local_vector.Import(gamma, importer, Insert);
		int curr_i = 0;
		if (has_north) {
			boundary_north = &local_vector[0][curr_i];
			curr_i += nx;
		}
		if (has_east) {
			boundary_east = &local_vector[0][curr_i];
			curr_i += ny;
		}
		if (has_south) {
			boundary_south = &local_vector[0][curr_i];
			curr_i += nx;
		}
		if (has_west) {
			boundary_west = &local_vector[0][curr_i];
		}

        // solve
        solve();

		local_vector.PutScalar(0.0);
		curr_i = 0;
		if (has_north) {
            for(int i = 0; i < nx; i++){
                local_vector[0][curr_i] = (*u)[0][i];
				curr_i++;
			}
        }
		if (has_east) {
            for(int i = 0; i < ny; i++){
                local_vector[0][curr_i] = (*u)[0][(i+1)*nx-1];
				curr_i++;
			}
		}
		if (has_south) {
            for(int i = 0; i < nx; i++){
                local_vector[0][curr_i] = (*u)[0][nx*(ny-1)+i];
				curr_i++;
			}
		}
		if (has_west) {
            for(int i = 0; i < ny; i++){
                local_vector[0][curr_i] = (*u)[0][i*nx];
				curr_i++;
			}
		}
        diff.Export(local_vector,importer,Add);
        diff.Update(-2,gamma,1);
	}
};

namespace Belos
{
class FuncWrap// : virtual Belos::Operator
{
    public:
	RCP<vector_type> b;
    Domain* d;
    FuncWrap(RCP<vector_type> b,Domain* d){
        this->b = b;
        this->d = d;
    }
    void Apply(const vector_type& x, vector_type& y){
        d->solveWithInterface(x,y);
        y.Update(-1,*b,1);
    }
    static bool HasApplyTranspose(const FuncWrap& wrapper){
        return false;
    }
    //static bool this_type_is_missing_a_specialization(){
    //    return false;
    //}
	static void Apply(const FuncWrap &wrapper, const vector_type &x, vector_type &y,
	                  Belos::ETrans trans = Belos::ETrans::NOTRANS)
	{
        Domain* d_ptr = wrapper.d;
        RCP<vector_type> b_ptr = wrapper.b;
        d_ptr->solveWithInterface(x,y);
        y.Update(-1,*b_ptr,1);
		//wrapper.Apply(x, y);
	}
};
template<>
void OperatorTraits<double, vector_type, FuncWrap>::Apply(const FuncWrap &wrapper,
                                                          const vector_type &x, vector_type &y,
                                                          Belos::ETrans trans)
{
	Domain *         d_ptr = wrapper.d;
	RCP<vector_type> b_ptr = wrapper.b;
	d_ptr->solveWithInterface(x, y);
	y.Update(-1, *b_ptr, 1);
	// wrapper.Apply(x, y);
};
}
matrix_type *generate2DLaplacian(const Teuchos::RCP<map_type> Map, const int nx, const int ny,
                                   const double h_x, const double h_y);

// the functions that we are using
double ffun(double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); }
double gfun(double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); }
// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
    using Belos::FuncWrap;
	using Teuchos::RCP;
	using Teuchos::rcp;

	MPI_Init(&argc, &argv);
	Epetra_MpiComm Comm(MPI_COMM_WORLD);

    int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);
    // comunicator for subdomain
    // for now, subdomains are going to be single threaded
	MPI_Comm subdomain_comm_raw;
	MPI_Comm_split(MPI_COMM_WORLD, my_global_rank, 0, &subdomain_comm_raw);
	Epetra_MpiComm subdomain_comm(subdomain_comm_raw);


	// parse input
	args::ArgumentParser  parser("");
	args::HelpFlag        help(parser, "help", "Display this help menu", {'h', "help"});
	args::Positional<int> d_x(parser, "d_x", "number of domains in the x direction");
	args::Positional<int> d_y(parser, "d_y", "number of domains in the y direction");
	args::Positional<int> n_x(parser, "n_x", "number of cells in the x direction, in each domain");
	args::Positional<int> n_y(parser, "n_y", "number of cells in the y direction, in each domain");
	args::ValueFlag<std::string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                                 {'m'});
	args::ValueFlag<std::string> f_s(parser, "solution filename",
	                                 "the file to write the solution to", {'s'});
	args::Flag f_cg(parser, "cg", "use conjugate gradient for solving gamma values", {"cg"});

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
	int    num_domains_x = args::get(d_x);
	int    num_domains_y = args::get(d_y);
	int    nx            = args::get(n_x);
	int    ny            = args::get(n_y);
	double h_x           = 1.0 / (nx * num_domains_x);
	double h_y           = 1.0 / (ny * num_domains_y);
	int    domain_x      = my_global_rank % num_domains_x;
	int    domain_y      = my_global_rank / num_domains_x;

	if (num_domains_x * num_domains_y != num_procs) {
		std::cerr << "number of domains must be equal to the number of processes\n";
        return 1;
    }

	// create a map and matrix
	RCP<map_type>    Map = rcp(new map_type(nx * ny, 0, subdomain_comm));
	RCP<matrix_type> A   = rcp(generate2DLaplacian(Map, nx, ny, h_x, h_y));

	// Generate RHS vector
	RCP<vector_type> f     = rcp(new vector_type(*Map, 1, false));
	RCP<vector_type> exact = rcp(new vector_type(*Map, 1, false));
	{
		// Use local indices to access the entries of f_data.
		const int localLength      = f->MyLength();
		int       NumMyElements    = Map->NumMyElements();
		int *     MyGlobalElements = 0;
		Map->MyGlobalElementsPtr(MyGlobalElements);
		// generate rhs
		for (int i = 0; i < NumMyElements; i++) {
			int    global_i = MyGlobalElements[i];
			int    index_x  = domain_x * nx + global_i % nx;
			int    index_y  = domain_y * ny + (global_i - index_x) / nx;
			double x        = h_x / 2.0 + 1.0 * index_x / (nx * num_domains_x);
			double y        = h_y / 2.0 + 1.0 * index_y / (ny * num_domains_y);
			(*f)[0][i]      = ffun(x, y);
			(*exact)[0][i]  = gfun(x, y);
			if (index_x == 0) {
				(*f)[0][i] += -2.0 / (h_x * h_x) * gfun(0.0, y);
			}
			if (index_x == num_domains_x * nx - 1) {
				(*f)[0][i] += -2.0 / (h_x * h_x) * gfun(1.0, y);
			}
			if (index_y == 0) {
				(*f)[0][i] += -2.0 / (h_y * h_y) * gfun(x, 0.0);
			}
			if (index_y == num_domains_y * ny - 1) {
				(*f)[0][i] += -2.0 / (h_y * h_y) * gfun(x, 1.0);
			}
		}
	}

	// create a domain
	Domain d(A, f, nx, ny, h_x, h_y);
    // generate indices for gamma vector
	int num_interface_points = 0;
    // north
	if (domain_y != num_domains_y - 1) {
		num_interface_points += nx;
		d.has_north = true;
	}
    // east
	if (domain_x != num_domains_x - 1) {
		num_interface_points += ny;
		d.has_east = true;
	}
    // south
	if (domain_y != 0) {
		num_interface_points += nx;
		d.has_south = true;
	}
    // west
	if (domain_x != 0) {
		num_interface_points += ny;
		d.has_west = true;
	}
	std::vector<int> global_i(num_interface_points);

	RCP<vector_type> u          = d.u;
	int              ns_start_i = num_domains_y * (num_domains_x - 1) * ny;
	int              curr_i     = 0;
	// north
	if (d.has_north) {
		int curr_global_i = (domain_y * num_domains_x + domain_x) * nx + ns_start_i;
		for (int i = 0; i < nx; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	// east
	if (d.has_east) {
		int curr_global_i = (domain_x * num_domains_y + domain_y) * ny;
		for (int i = 0; i < ny; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	// south
	if (d.has_south) {
		int curr_global_i = ((domain_y - 1) * num_domains_x + domain_x) * nx + ns_start_i;
		for (int i = 0; i < nx; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	// west
	if (d.has_west) {
		int curr_global_i = ((domain_x - 1) * num_domains_y + domain_y) * ny;
		for (int i = 0; i < ny; i++) {
			global_i[curr_i] = curr_global_i;
			curr_global_i++;
			curr_i++;
		}
	}
	// create the map
	int num_global_elements = nx * num_domains_x * (num_domains_y - 1) + ny * num_domains_y * (num_domains_x - 1);

	RCP<map_type> diff_map = rcp(new map_type(num_global_elements, 0, Comm));

	if (num_domains_x * num_domains_y == 1) {
		d.domain_map = rcp(new map_type(1, 0, Comm));
	} else {
		d.domain_map = rcp(new map_type(-1, num_interface_points, &global_i[0], 0, Comm));
	}

	RCP<vector_type> gamma    = rcp(new vector_type(*diff_map, 1));
	RCP<vector_type> diff     = rcp(new vector_type(*diff_map, 1));
	if (num_domains_x * num_domains_y != 1) {
        //do iterative solve
		RCP<vector_type> b = rcp(new vector_type(*diff_map, 1));
		d.solveWithInterface(*gamma, *b);
        b->Print(std::cout);
        sleep(1);
		RCP<FuncWrap> wrapper = rcp(new FuncWrap(b, &d));
		Belos::LinearProblem<double, vector_type, FuncWrap> problem(wrapper, gamma, b);
        std::cout<< problem.setProblem()<<"\n";
        Teuchos::ParameterList belosList;
        belosList.set("Block Size",1);
        belosList.set("Maximum Iterations",1000);
        belosList.set("Convergence Tolerance",10e-7);
		int verbosity = Belos::Errors + Belos::StatusTestDetails + Belos::Warnings
		                + Belos::TimingDetails + Belos::Debug;
        belosList.set("Verbosity",verbosity);
        Belos::OutputManager<double> my_om();
		RCP < Belos::SolverManager<double, vector_type, FuncWrap>> solver
		= rcp(new Belos::BlockCGSolMgr<double, vector_type, FuncWrap>(rcp(&problem, false),
		                                                              rcp(&belosList, false)));
        solver->solve();
        
	}

	//d.domain_map->Print(std::cout);
    //diff_map->Print(std::cout);
	d.solveWithInterface(*gamma,*diff);
    map_type err_map(-1,1,0,Comm);
    Epetra_Vector exact_norm(err_map);
	Epetra_Vector diff_norm(err_map);
	exact->Norm2(&exact_norm[0]);
	{
		// Use local indices to access the entries of f_data.
		const int localLength      = f->MyLength();
		int       NumMyElements    = Map->NumMyElements();
		int *     MyGlobalElements = 0;
		Map->MyGlobalElementsPtr(MyGlobalElements);
		for (int i = 0; i < NumMyElements; i++) {
			(*exact)[0][i] -= (*u)[0][i];
		}
	}
	exact->Norm2(&diff_norm[0]);
	double global_diff_norm;
	double global_exact_norm;
	diff_norm.Norm2(&global_diff_norm);
	exact_norm.Norm2(&global_exact_norm);
	if (my_global_rank == 0) {
		std::cout << global_diff_norm / global_exact_norm << "\n";
	}

	MPI_Finalize();
	return 0;
}

matrix_type *generate2DLaplacian(const Teuchos::RCP<map_type> Map, const int nx,
                                      const int ny, const double h_x, const double h_y)
{
	using Teuchos::RCP;
	matrix_type *Matrix = new matrix_type(Copy, *Map, 5);

	int  NumMyElements    = Map->NumMyElements();
	int *MyGlobalElements = 0;
	Map->MyGlobalElementsPtr(MyGlobalElements);

	int                 left, right, lower, upper;
	std::vector<double> Values(4);
	std::vector<int>    Indices(4);

	//    e
	//  b a c
	//    d
	for (int i = 0; i < NumMyElements; ++i) {
		int NumEntries = 0;
		// determine index of of neighbors in row
		int global_i = MyGlobalElements[i];
		int ix, iy;
		ix = global_i % nx;
		iy = (global_i - ix) / nx;

		if (ix == 0)
			left = -1;
		else
			left = global_i - 1;
		if (ix == nx - 1)
			right = -1;
		else
			right = global_i + 1;
		if (iy == 0)
			lower = -1;
		else
			lower = global_i - nx;
		if (iy == ny - 1)
			upper = -1;
		else
			upper = global_i + nx;

		double diag = -2.0 / (h_y * h_y) - 2.0 / (h_x * h_x);

		if (left != -1) {
			Indices[NumEntries] = left;
			Values[NumEntries]  = 1.0 / (h_y * h_y);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_y * h_y);
		}
		if (right != -1) {
			Indices[NumEntries] = right;
			Values[NumEntries]  = 1.0 / (h_y * h_y);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_y * h_y);
		}
		if (lower != -1) {
			Indices[NumEntries] = lower;
			Values[NumEntries]  = 1.0 / (h_x * h_x);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_x * h_x);
		}
		if (upper != -1) {
			Indices[NumEntries] = upper;
			Values[NumEntries]  = 1.0 / (h_x * h_x);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_x * h_x);
		}
		// put the off-diagonal entries
		Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);

		// Put in the diagonal entry
		Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements + i);
	}
	Matrix->FillComplete();
	Matrix->OptimizeStorage();

	return (Matrix);
}

