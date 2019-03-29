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

#ifndef GMGMatOp_H
#define GMGMatOp_H
#include <GMG/Operator.h>
#include <petscmat.h>
#include <memory>
namespace GMG
{
/**
 * @brief Wrapper for PETSc Matrix
 */
class MatOp : public Operator
{
	private:
	/**
	 * @brief PETSc Matrix object
	 */
	PW<Mat> matrix;

	public:
	/**
	 * @brief Crate new MatOp
	 *
	 * @param matrix the PETSc matrix
	 */
	MatOp(PW<Mat> matrix)
	{
		this->matrix = matrix;
	}
	/**
	 * @brief Perform matrix/vector multiply.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(PW<Vec> x, PW<Vec> b) const
	{
		MatMult(matrix, x, b);
	}
};
} // namespace GMG
#endif
