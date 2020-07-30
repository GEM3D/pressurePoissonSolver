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

#ifndef PETSCMATOP_H
#define PETSCMATOP_H
#include <Thunderegg/Operators/Operator.h>
#include <Thunderegg/SchurHelper.h>
/**
 * @brief Base class for operators
 */
template <size_t D> class PetscMatOp : public Operator<D>
{
	private:
	PW<Mat> A;

	public:
	PetscMatOp(PW<Mat> A)
	{
		this->A = A;
	}
	/**
	 * @brief Apply Petsc matrix
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
	{
		const PetscVector<D> *x_vec = dynamic_cast<const PetscVector<D> *>(x.get());
		PetscVector<D> *      b_vec = dynamic_cast<PetscVector<D> *>(b.get());
		if (x_vec == nullptr || b_vec == nullptr) { throw 3; }
		MatMult(A, x_vec->vec, b_vec->vec);
	}
};
#endif
