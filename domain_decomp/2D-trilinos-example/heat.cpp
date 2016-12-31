#include <Galeri_CrsMatrices.h>
#include <Galeri_Maps.h>
#include <Galeri_Utils.h>
#include <Epetra_MpiComm.h>
#include <mpi.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

using namespace Galeri;
	using Teuchos::RCP;

Epetra_CrsMatrix *Cross2D(const RCP<Epetra_Map> Map, const int nx, const int ny, const double a,
                          const double b, const double c, const double d, const double e);
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
	GaleriList.set("nx", nx);
	GaleriList.set("ny", ny);
	GaleriList.set("lx", 1.0);
	GaleriList.set("ly", 1.0);

	// Create the map and matrix using the parameter list for a 2D Laplacian.
	RCP<Epetra_Map>       Map    = rcp(CreateMap("Cartesian2D", Comm, GaleriList));
	RCP<Epetra_CrsMatrix> Matrix = rcp(Cross2D(Map,nx,ny,4.0,-1.0,-1.0,-1.0,-1.0));

    // Create grid for rhs
    
	// Print out the map and matrices
	Map->Print(std::cout);
	Matrix->Print(std::cout);

	MPI_Finalize();
	return 0;
}

void GetNeighboursCartesian2d(const int i, const int nx, const int ny, int &left, int &right,
                              int &lower, int &upper, int &left2, int &right2, int &lower2,
                              int &upper2)
{
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

	if (ix <= 1)
		left2 = -1;
	else
		left2 = i - 2;
	if (ix >= nx - 2)
		right2 = -1;
	else
		right2 = i + 2;
	if (iy <= 1)
		lower2 = -1;
	else
		lower2 = i - 2 * nx;
	if (iy >= ny - 2)
		upper2 = -1;
	else
		upper2 = i + 2 * nx;
}

Epetra_CrsMatrix *Cross2D(const RCP<Epetra_Map> Map, const int nx, const int ny, const double a,
                          const double b, const double c, const double d, const double e)
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
		GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, left, right, lower, upper);

		double diag = a;

		if (left != -1) {
			Indices[NumEntries] = left;
			Values[NumEntries]  = b;
			++NumEntries;
		} else {
			diag += 1.0;
		}
		if (right != -1) {
			Indices[NumEntries] = right;
			Values[NumEntries]  = c;
			++NumEntries;
		} else {
			diag += 1.0;
		}
		if (lower != -1) {
			Indices[NumEntries] = lower;
			Values[NumEntries]  = d;
			++NumEntries;
		} else {
			diag += 1.0;
		}
		if (upper != -1) {
			Indices[NumEntries] = upper;
			Values[NumEntries]  = e;
			++NumEntries;
		} else {
			diag += 1.0;
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

