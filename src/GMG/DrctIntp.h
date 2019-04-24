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

#ifndef GMGDrctIntp_H
#define GMGDrctIntp_H
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
template <size_t D> class DrctIntp : public Interpolator<D>
{
	private:
	/**
	 * @brief The coarser set of domains
	 */
	std::shared_ptr<DomainCollection<D>> coarse_dc;
	/**
	 * @brief The finer set of domains
	 */
	std::shared_ptr<DomainCollection<D>> fine_dc;
	/**
	 * @brief The comm package between the levels.
	 */
	std::shared_ptr<InterLevelComm<D>> ilc;

	public:
	/**
	 * @brief Create new DrctIntp object.
	 *
	 * @param coarse_dc the coarser set of domains.
	 * @param fine_dc the finer set of domains.
	 * @param ilc the comm package between the levels.
	 */
	DrctIntp(std::shared_ptr<DomainCollection<D>> coarse_dc,
	         std::shared_ptr<DomainCollection<D>> fine_dc, std::shared_ptr<InterLevelComm<D>> ilc);
	/**
	 * @brief Interpolate from the finer level to the coarser level.
	 *
	 * @param coarse the input vector from the coarser level
	 * @param fine the output vector for the finer level
	 */
	void interpolate(std::shared_ptr<const Vector<D>> coarse, std::shared_ptr<Vector<D>> fine)const;
};
template <size_t D>
inline DrctIntp<D>::DrctIntp(std::shared_ptr<DomainCollection<D>> coarse_dc,
                             std::shared_ptr<DomainCollection<D>> fine_dc,
                             std::shared_ptr<InterLevelComm<D>>   ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}

template <size_t D>
inline void DrctIntp<D>::interpolate(std::shared_ptr<const Vector<D>> coarse,
                                     std::shared_ptr<Vector<D>>       fine) const
{
	// scatter
	std::shared_ptr<Vector<D>> coarse_local = ilc->getNewCoarseDistVec();
	ilc->scatter(coarse_local, coarse);

	for (auto data : ilc->getFineDomains()) {
		Domain<D> &d          = *data.d;
		LocalData<D> coarse_local_data = coarse_local->getLocalData(data.local_index);
		LocalData<D> fine_data         = fine->getLocalData(d.id_local);

		if (d.hasCoarseParent()) {
			Orthant<D>         orth = d.oct_on_parent;
			std::array<int, D> starts;
			for (size_t i = 0; i < D; i++) {
				starts[i]  = orth.isOnSide(2 * i) ? 0 : coarse_local_data.getLengths()[i];
			}

			nested_loop<D>(fine_data.getStart(), fine_data.getEnd(),
			               [&](const std::array<int, D> &coord) {
				               std::array<int, D> coarse_coord;
				               for (size_t x = 0; x < D; x++) {
					               coarse_coord[x] = (coord[x] + starts[x]) / 2;
				               }
				               fine_data[coord] += coarse_local_data[coarse_coord];
			               });
		} else {
			nested_loop<D>(
			fine_data.getStart(), fine_data.getEnd(),
			[&](const std::array<int, D> &coord) { fine_data[coord] += coarse_local_data[coord]; });
		}
	}
}
} // namespace GMG
#endif
