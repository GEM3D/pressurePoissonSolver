#ifndef GENERATELAPLACIAN_H
#define GENERATELAPLACIAN_H
#include "MyTypeDefs.h"
matrix_type *generate2DLaplacian(const Teuchos::RCP<map_type> Map, const int nx, const int ny,
                                 const double h_x, const double h_y)
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
			Values[NumEntries]  = 1.0 / (h_x * h_x);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_x * h_x);
		}
		if (right != -1) {
			Indices[NumEntries] = right;
			Values[NumEntries]  = 1.0 / (h_x * h_x);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_x * h_x);
		}
		if (lower != -1) {
			Indices[NumEntries] = lower;
			Values[NumEntries]  = 1.0 / (h_y * h_y);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_y * h_y);
		}
		if (upper != -1) {
			Indices[NumEntries] = upper;
			Values[NumEntries]  = 1.0 / (h_y * h_y);
			++NumEntries;
		} else {
			diag -= 1.0 / (h_y * h_y);
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
#endif
