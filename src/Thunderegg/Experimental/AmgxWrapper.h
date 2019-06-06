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

#ifndef AMGXWRAPPER_H
#define AMGXWRAPPER_H
#include "DomainCollection.h"
#include "amgx_c.h"
#include <petscmat.h>
#include <iostream>
#include <string>
class AmgxWrapper
{
	private:
	static void print_callback(const char *msg, int length)
	{
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0) std::cout << msg;
	}
	// library handles
	MPI_Comm              AMGX_MPI_COMM;
	AMGX_config_handle    cfg;
	AMGX_resources_handle rsrc;
	AMGX_matrix_handle    gA;
	AMGX_vector_handle    gb;
	AMGX_vector_handle    gx;
	AMGX_solver_handle    solver;
	int               num_rows;
        int nrings = 0;

	public:
	AmgxWrapper(std::string filename);
	void setMatrix(Mat A);
	~AmgxWrapper();
	void solve(Vec x, Vec b);
};
#endif
