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

#include "IfaceMatrixHelper.h"
#include "IMStencilHelper.h"
#include <iostream>
using namespace std;
IfaceMatrixHelper::IfaceMatrixHelper(DomainCollection dc) { this->dc = dc; }
PW_explicit<Mat> IfaceMatrixHelper::formCRSMatrix(double lambda)
{
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int n           = dc.n;
	int local_size  = dc.domains.size() * dc.n * dc.n + dc.ifaces.size() * dc.n;
	int global_size = dc.domains.size() * dc.n * dc.n + dc.ifaces.size() * dc.n;
	int iface_start = dc.domains.size() * dc.n * dc.n;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19, nullptr, 19, nullptr);

	for (auto &p : dc.domains) {
		Domain &d     = p.second;
		int     n     = d.n;
		double  h_x   = d.x_length / n;
		double  h_y   = d.y_length / n;
		int     start = n * n * d.id_global;

		// center coeffs
		double coeff = -2.0 / (h_x * h_x) - 2.0 / (h_y * h_y) + lambda;
		for (int y_i = 0; y_i < n; y_i++) {
			for (int x_i = 0; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				MatSetValues(A, 1, &row, 1, &row, &coeff, ADD_VALUES);
			}
		}
		// north coeffs
		coeff = 1.0 / (h_y * h_y);
		for (int y_i = 0; y_i < n - 1; y_i++) {
			for (int x_i = 0; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i + n * (y_i + 1);
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		// south coeffs
		for (int y_i = 1; y_i < n; y_i++) {
			for (int x_i = 0; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i + n * (y_i - 1);
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		coeff = 1.0 / (h_x * h_x);
		// east coeffs
		for (int y_i = 0; y_i < n; y_i++) {
			for (int x_i = 0; x_i < n - 1; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i + 1 + n * y_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		// west coeffs
		for (int y_i = 0; y_i < n; y_i++) {
			for (int x_i = 1; x_i < n; x_i++) {
				int row = start + x_i + n * y_i;
				int col = start + x_i - 1 + n * y_i;
				MatSetValues(A, 1, &row, 1, &col, &coeff, ADD_VALUES);
			}
		}
		// boundaries
		Side s = Side::north;
		do {
			IMStencilHelper *sh = getIMStencilHelper(d, s, iface_start);
			for (int i = 0; i < n; i++) {
				int     row    = sh->row(i);
				int     size   = sh->size(i);
				double *coeffs = sh->coeffs(i);
				int *   cols   = sh->cols(i);
				MatSetValues(A, 1, &row, size, cols, coeffs, ADD_VALUES);
			}
			delete sh;
			s++;
		} while (s != Side::north);
	}
	for (auto &p : dc.ifaces) {
		Iface &          iface = p.second;
		IMStencilHelper *sh    = getIfaceStencilHelper(iface, iface_start, n);
		for (int i = 0; i < n; i++) {
			int     row    = sh->row(i);
			int     size   = sh->size(i);
			double *coeffs = sh->coeffs(i);
			int *   cols   = sh->cols(i);
			MatSetValues(A, 1, &row, size, cols, coeffs, ADD_VALUES);
		}
		delete sh;
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
