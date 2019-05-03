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

#ifndef SCHURWRAPOP_H
#define SCHURWRAPOP_H
#include <Thunderegg/Operators/Operator.h>
#include <Thunderegg/SchurHelper.h>
/**
 * @brief Base class for operators
 */
template <size_t D> class SchurWrapOp : public Operator<D - 1>
{
	private:
	std::shared_ptr<DomainCollection<D>> dc;
	std::shared_ptr<SchurHelper<D>>      sh;

	public:
	SchurWrapOp(std::shared_ptr<DomainCollection<D>> dc, std::shared_ptr<SchurHelper<D>> sh)
	{
		this->dc = dc;
		this->sh = sh;
	}
	/**
	 * @brief Apply Schur matrix
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D - 1>> x, std::shared_ptr<Vector<D - 1>> b) const
	{
		auto f = dc->getNewDomainVec();
		auto u = dc->getNewDomainVec();
		sh->solveAndInterpolateWithInterface(f, u, x, b);
	}
};
#endif
