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
#include <DomainCollection.h>
#include <GMG/InterLevelComm.h>
#include <GMG/Interpolator.h>
#include <memory>
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse vector to fine vector. (This is
 * O(0) accuracy, need to replace with bi-linear interpolation.
 */
class TriLinIntp : public Interpolator
{
	private:
	/**
	 * @brief The coarser set of domains
	 */
	std::shared_ptr<DomainCollection<3>> coarse_dc;
	/**
	 * @brief The finer set of domains
	 */
	std::shared_ptr<DomainCollection<3>> fine_dc;
	/**
	 * @brief The comm package between the levels.
	 */
	std::shared_ptr<InterLevelComm<3>> ilc;

	public:
	/**
	 * @brief Create new TriLinIntp object.
	 *
	 * @param coarse_dc the coarser set of domains.
	 * @param fine_dc the finer set of domains.
	 * @param ilc the comm package between the levels.
	 */
	TriLinIntp(std::shared_ptr<DomainCollection<3>> coarse_dc,
	           std::shared_ptr<DomainCollection<3>> fine_dc,
	           std::shared_ptr<InterLevelComm<3>>   ilc);
	/**
	 * @brief Interpolate from the finer level to the coarser level.
	 *
	 * @param coarse the input vector from the coarser level
	 * @param fine the output vector for the finer level
	 */
	void interpolate(PW<Vec> coarse, PW<Vec> fine) const;
};
} // namespace GMG
#endif
