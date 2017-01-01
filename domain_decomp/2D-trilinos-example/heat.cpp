#include <Galeri_CrsMatrices.h>
#include <Galeri_Maps.h>
#include <Galeri_Utils.h>
#include <Epetra_MpiComm.h>
#include <mpi.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

using namespace Galeri;
using Teuchos::RCP;

Epetra_CrsMatrix *generate2DLaplacian(const RCP<Epetra_Map> Map, const int nx, const int ny, const double h_x, const double h_y);

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
    int nx = 5;
    int ny = 10;
	double h_x = 1.0 / nx;
	double h_y = 1.0 / ny;
	GaleriList.set("nx", nx);
	GaleriList.set("ny", ny);
	GaleriList.set("lx", 1.0);
	GaleriList.set("ly", 1.0);

	// Create the map and matrix using the parameter list for a 2D Laplacian.
	RCP<Epetra_Map>       Map    = rcp(CreateMap("Cartesian2D", Comm, GaleriList));

	RCP<Epetra_CrsMatrix> Matrix = rcp(generate2DLaplacian(Map,nx,ny,h_x,h_y));

    // Generate RHS vector
	Epetra_Vector f(*Map,nx*ny);
	{
		// Use local indices to access the entries of f_data.
		const int localLength      = f.MyLength();
		int       NumMyElements    = Map->NumMyElements();
		int *     MyGlobalElements = 0;
		Map->MyGlobalElementsPtr(MyGlobalElements);
		cout << MyGlobalElements << "\n";
		for (int i = 0; i < NumMyElements; i++) {
			int    global_i = MyGlobalElements[i];
			int    index_x  = i % nx;
			int    index_y  = (i - index_x) / nx;
			double x        = h_x / 2.0 + 1.0 * index_x / nx;
			double y        = h_y / 2.0 + 1.0 * index_y / ny;
			f[i]            = ffun(x, y);
		}
	}

	// Print out the map and matrices
	Map->Print(std::cout);
	Matrix->Print(std::cout);
	f.Print(std::cout);

	MPI_Finalize();
	return 0;
}

Epetra_CrsMatrix *generate2DLaplacian(const RCP<Epetra_Map> Map, const int nx, const int ny, const double h_x, const double h_y)
{
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
		int ix, iy;
		ix = i % nx;
		iy = (i - ix) / nx;

		if (ix == 0)
			left = -1;
		else
			left = i - 1;
		if (ix == nx - 1)
			right = -1;
		else
			right = i + 1;
		if (iy == 0)
			lower = -1;
		else
			lower = i - nx;
		if (iy == ny - 1)
			upper = -1;
		else
			upper = i + nx;


		double diag = -2.0/(h_y*h_y)-2.0/(h_x*h_x);

		if (left != -1) {
			Indices[NumEntries] = left;
			Values[NumEntries]  = 1.0/(h_y*h_y);
			++NumEntries;
		} else {
			diag -= 1.0/(h_y*h_y);
		}
		if (right != -1) {
			Indices[NumEntries] = right;
			Values[NumEntries]  = 1.0/(h_y*h_y);
			++NumEntries;
		} else {
			diag -= 1.0/(h_y*h_y);
		}
		if (lower != -1) {
			Indices[NumEntries] = lower;
			Values[NumEntries]  = 1.0/(h_x*h_x);
			++NumEntries;
		} else {
			diag -= 1.0/(h_x*h_x);
		}
		if (upper != -1) {
			Indices[NumEntries] = upper;
			Values[NumEntries]  = 1.0/(h_x*h_x);
			++NumEntries;
		} else {
			diag -= 1.0/(h_x*h_x);
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

