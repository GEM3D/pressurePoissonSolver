#include <Galeri_CrsMatrices.h>
#include <Galeri_Maps.h>
#include <Galeri_Utils.h>
#include <Epetra_MpiComm.h>
#include <mpi.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>

using namespace Galeri;

Epetra_CrsMatrix *generate2DLaplacian(const Teuchos::RCP<Epetra_Map> Map, const int nx,
                                      const int ny, const double h_x, const double h_y);

// the functions that we are using
double ffun(double x, double y) { return -5 * M_PI * M_PI * sin(M_PI * x) * cos(2 * M_PI * y); }
double gfun(double x, double y) { return sin(M_PI * x) * cos(2 * M_PI * y); }
// =========== //
// main driver //
// =========== //

int main(int argv, char *argc[])
{
	using Teuchos::RCP;
	using Teuchos::rcp;

	MPI_Init(&argv, &argc);
	Epetra_MpiComm Comm(MPI_COMM_WORLD);

	// Create a parameter list
	Teuchos::ParameterList GaleriList;

	// Set the number of discretization points in the x and y direction.
	int    nx  = 256;
	int    ny  = 256;
	double h_x = 1.0 / nx;
	double h_y = 1.0 / ny;
	GaleriList.set("nx", nx);
	GaleriList.set("ny", ny);
	GaleriList.set("lx", 1.0);
	GaleriList.set("ly", 1.0);

	// Create the map and matrix using the parameter list for a 2D Laplacian.
	RCP<Epetra_Map> Map = rcp(new Epetra_Map(nx * ny, 0, Comm));

	RCP<Epetra_CrsMatrix> A = rcp(generate2DLaplacian(Map, nx, ny, h_x, h_y));

	// Generate RHS vector
	RCP<Epetra_MultiVector> f     = rcp(new Epetra_MultiVector(*Map, 1));
	RCP<Epetra_MultiVector> u     = rcp(new Epetra_MultiVector(*Map, 1));
	RCP<Epetra_MultiVector> exact = rcp(new Epetra_MultiVector(*Map, 1));
	{
		// Use local indices to access the entries of f_data.
		const int localLength      = f->MyLength();
		int       NumMyElements    = Map->NumMyElements();
		int *     MyGlobalElements = 0;
		Map->MyGlobalElementsPtr(MyGlobalElements);
		cout << MyGlobalElements << "\n";
		for (int i = 0; i < NumMyElements; i++) {
			int    global_i = MyGlobalElements[i];
			int    index_x  = global_i % nx;
			int    index_y  = (global_i - index_x) / nx;
			double x        = h_x / 2.0 + 1.0 * index_x / nx;
			double y        = h_y / 2.0 + 1.0 * index_y / ny;
			(*f)[0][i]      = ffun(x, y);
			(*exact)[0][i]  = gfun(x, y);
		}
		// add boundaries to rhs
		for (int i = 0; i < NumMyElements; i++) {
			int    global_i = MyGlobalElements[i];
			int    index_x  = global_i % nx;
			int    index_y  = (global_i - index_x) / nx;
			double x        = h_x / 2.0 + 1.0 * index_x / nx;
			double y        = h_y / 2.0 + 1.0 * index_y / ny;
			if (index_x == 0) {
				(*f)[0][i] += -2.0 / (h_x * h_x) * gfun(0.0, y);
			}
			if (index_x == nx - 1) {
				(*f)[0][i] += -2.0 / (h_x * h_x) * gfun(1.0, y);
			}
			if (index_y == 0) {
				(*f)[0][i] += -2.0 / (h_y * h_y) * gfun(x, 0.0);
			}
			if (index_y == ny - 1) {
				(*f)[0][i] += -2.0 / (h_y * h_y) * gfun(x, 1.0);
			}
		}
	}

	// Print out the map and matrices
	// Map->Print(std::cout);
	// A->Print(std::cout);
	// f->Print(std::cout);
	// exact->Print(std::cout);

	// Create solver interface with Amesos2 factory method
	RCP<Amesos2::Solver<Epetra_CrsMatrix, Epetra_MultiVector> > solver
	= Amesos2::create<Epetra_CrsMatrix, Epetra_MultiVector>("KLU2", A, u, f);
	solver->solve();

	double exact_norm;
	exact->Norm2(&exact_norm);
	{
		// Use local indices to access the entries of f_data.
		const int localLength      = f->MyLength();
		int       NumMyElements    = Map->NumMyElements();
		int *     MyGlobalElements = 0;
		Map->MyGlobalElementsPtr(MyGlobalElements);
		cout << MyGlobalElements << "\n";
		for (int i = 0; i < NumMyElements; i++) {
			(*exact)[0][i] -= (*u)[0][i];
		}
	}
	// u->Print(std::cout);
	double diff_norm;
	exact->Norm2(&diff_norm);
	cout << diff_norm / exact_norm << "\n";

	MPI_Finalize();
	return 0;
}

Epetra_CrsMatrix *generate2DLaplacian(const Teuchos::RCP<Epetra_Map> Map, const int nx,
                                      const int ny, const double h_x, const double h_y)
{
	using Teuchos::RCP;
	Epetra_CrsMatrix *Matrix = new Epetra_CrsMatrix(Copy, *Map, 5);

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

