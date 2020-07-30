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

#ifndef GMGFFTBlockJacobiSmoother_H
#define GMGFFTBlockJacobiSmoother_H
#include <Thunderegg/SchurHelper.h>
namespace GMG
{
/**
 * @brief A block Jacobi smoother that uses FFTW solves on each patch. Implemented using the
 * SchurHelper class.
 */
template <size_t D> class FFTBlockJacobiSmoother : public Smoother<D>
{
	private:
	/**
	 * @brief point to the SchurHelper object.
	 */
	std::shared_ptr<SchurHelper<D>> sh;

	public:
	/**
	 * @brief Create new smoother with SchurHelper object
	 *
	 * @param sh pointer to the SchurHelper object
	 */
	FFTBlockJacobiSmoother(std::shared_ptr<SchurHelper<D>> sh)
	{
		this->sh = sh;
	}
	/**
	 * @brief Run an iteration of smoothing.
	 *
	 * @param f the RHS vector
	 * @param u the solution vector, updated upon return.
	 */
	void smooth(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u) const
	{
		sh->solveWithSolution(f, u);
	}
};
} // namespace GMG
#endif
