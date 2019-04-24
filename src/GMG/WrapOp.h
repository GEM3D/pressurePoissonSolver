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

#ifndef GMGWrapOp_H
#define GMGWrapOp_H
#include <GMG/Operator.h>
#include <memory>
namespace GMG
{
/**
 * @brief Wrapper for Matrix free operation
 */
template <size_t D> class WrapOp : public Operator<D>
{
	private:
	/**
	 * @brief PETSc Matrix object
	 */
	std::shared_ptr<SchurHelper<D>> helper;

	public:
	/**
	 * @brief Crate new WrapOp
	 *
	 * @param matrix the PETSc matrix
	 */
	WrapOp(std::shared_ptr<SchurHelper<D>> helper)
	{
		this->helper = helper;
	}
	/**
	 * @brief Perform matrix/vector multiply.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
	{
		helper->apply(x, b);
	}
};
} // namespace GMG
#endif
