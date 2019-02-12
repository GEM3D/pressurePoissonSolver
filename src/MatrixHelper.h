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

#ifndef MATRIXHELPER_H
#define MATRIXHELPER_H
#include "DomainCollection.h"
#include <petscmat.h>
class MatrixHelper
{
	private:
	DomainCollection<3> dc;

	public:
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param dc the DomainCollection
	 * @param n number of cells in each direction for each domain
	 * @param h_x the x spacing
	 * @param h_y the y spacing
	 * @param comm the teuchos communicator
	 */
	MatrixHelper(DomainCollection<3> dc);

	/**
	 * @brief Form the matrix
	 *
	 * @return the formed matrix
	 */
	PW_explicit<Mat> formCRSMatrix(double lambda = 0);
};
#endif
