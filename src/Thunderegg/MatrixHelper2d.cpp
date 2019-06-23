/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#include "MatrixHelper2d.h"
#include "StencilHelper2d.h"
#include <iostream>
using namespace std;
MatrixHelper2d::MatrixHelper2d(std::shared_ptr<Domain<2>> domain)
{
	this->domain = domain;
}
PW_explicit<Mat> MatrixHelper2d::formCRSMatrix(double lambda)
{
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int nx          = domain->getNs()[0];
	int ny          = domain->getNs()[1];
	int local_size  = domain->getNumLocalPatches() * nx * ny;
	int global_size = domain->getNumGlobalPatches() * nx * ny;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19, nullptr, 19, nullptr);

	for (auto &pinfo : domain->getPatchInfoVector()) {
		double h_x   = pinfo->spacings[0];
		double h_y   = pinfo->spacings[1];
		int    start = nx * ny * pinfo->local_index;

		// center coeffs
		double coeff = -2.0 / (h_x * h_x) - 2.0 / (h_y * h_y) + lambda;
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i;
				MatSetValues(A, 1, &row, 1, &row, &coeff, ADD_VALUES);
			}
		}
		// north coeffs
		coeff = 1.0 / (h_y * h_y);
		for (int y_i = 0; y_i < ny - 1; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i;
				int col = start + x_i + nx * (y_i + 1);
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		// south coeffs
		for (int y_i = 1; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i;
				int col = start + x_i + nx * (y_i - 1);
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		coeff = 1.0 / (h_x * h_x);
		// east coeffs
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx - 1; x_i++) {
				int row = start + x_i + nx * y_i;
				int col = start + x_i + 1 + nx * y_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		// west coeffs
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 1; x_i < nx; x_i++) {
				int row = start + x_i + nx * y_i;
				int col = start + x_i - 1 + nx * y_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		// boundaries
		for (Side<2> s : Side<2>::getValues()) {
			StencilHelper *sh = getStencilHelper(*pinfo, s);
			for (int i = 0; i < sh->n; i++) {
				int     row    = sh->row(i);
				int     size   = sh->size(i);
				double *coeffs = sh->coeffs(i);
				int *   cols   = sh->cols(i);
				MatSetValues(A, 1, &row, size, cols, coeffs, ADD_VALUES);
			}
			delete sh;
		}
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	/*
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (domain->neumann&&rank==0) {
	    int           ncols;
	    const int *   cols;
	    const double *vals;
	    int row = 0;
	    MatGetRow(A, 0, &ncols, &cols, &vals);
	    vector<double> zeros(ncols);
	    zeros[0]=1;
	    MatSetValues(A, 1, &row, ncols, cols, &zeros[0], INSERT_VALUES);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	    MatRestoreRow(A, 0, &ncols, &cols, &vals);
	}
	*/
	return A;
}
