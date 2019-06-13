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

#ifndef GMGTriLinIntp_H
#define GMGTriLinIntp_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/Interpolator.h>
#include <memory>
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse vector to fine vector. (This is
 * O(0) accuracy, need to replace with bi-linear interpolation.
 */
class TriLinIntp : public Interpolator<3>
{
	private:
	/**
	 * @brief The coarser set of domains
	 */
	std::shared_ptr<Domain<3>> coarse_domain;
	/**
	 * @brief The finer set of domains
	 */
	std::shared_ptr<Domain<3>> fine_domain;
	/**
	 * @brief The comm package between the levels.
	 */
	std::shared_ptr<InterLevelComm<3>> ilc;

	public:
	/**
	 * @brief Create new TriLinIntp object.
	 *
	 * @param coarse_domain the coarser set of domains.
	 * @param fine_domain the finer set of domains.
	 * @param ilc the comm package between the levels.
	 */
	TriLinIntp(std::shared_ptr<Domain<3>> coarse_domain, std::shared_ptr<Domain<3>> fine_domain,
	           std::shared_ptr<InterLevelComm<3>> ilc);
	/**
	 * @brief Interpolate from the finer level to the coarser level.
	 *
	 * @param coarse the input vector from the coarser level
	 * @param fine the output vector for the finer level
	 */
	void interpolate(std::shared_ptr<const Vector<3>> coarse,
	                 std::shared_ptr<Vector<3>>       fine) const;
};
} // namespace GMG
#endif
