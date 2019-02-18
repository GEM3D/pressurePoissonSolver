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
#include "DomainCollection.h"
#include "InterLevelComm.h"
#include "Interpolator.h"
#include <memory>
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse vector to fine vector. (This is
 * O(0) accuracy, need to replace with bi-linear interpolation.
 */
template <size_t D> class DrctIntp : public Interpolator
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
	/**
	 * @brief pre-calculated powers of n
	 */
	std::array<int, D + 1> npow;

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
	void interpolate(PW<Vec> coarse, PW<Vec> fine) const;
};
template <size_t D>
inline DrctIntp<D>::DrctIntp(std::shared_ptr<DomainCollection<D>> coarse_dc,
                             std::shared_ptr<DomainCollection<D>> fine_dc,
                             std::shared_ptr<InterLevelComm<D>>   ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
	for (size_t i = 0; i <= D; i++) {
		npow[i] = (int) std::pow(fine_dc->getN(), i);
	}
}
template <size_t D> inline void DrctIntp<D>::interpolate(PW<Vec> coarse, PW<Vec> fine) const
{
	// get vectors
	double *u_fine;
	double *u_coarse;
	PW<Vec> coarse_tmp = ilc->getNewCoarseDistVec();
	// scatter
	PW<VecScatter> scatter = ilc->getScatter();
	VecScatterBegin(scatter, coarse, coarse_tmp, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, coarse, coarse_tmp, INSERT_VALUES, SCATTER_FORWARD);

	VecGetArray(fine, &u_fine);
	VecGetArray(coarse_tmp, &u_coarse);
	for (auto p : ilc->getFineDomains()) {
		Domain<D> &d          = *p.d;
		int        n          = d.n;
		int        coarse_idx = p.local_index * npow[D];
		int        fine_idx   = d.id_local * npow[D];

		Orthant<D> orth = d.oct_on_parent;
		if (d.id != d.parent_id) {
			std::array<int, D> strides;
			std::array<int, D> starts;
			for (size_t i = 0; i < D; i++) {
				strides[i] = npow[i];
				starts[i]  = orth.isOnSide(2 * i) ? 0 : n;
			}
			for (int i = 0; i < npow[D]; i++) {
				std::array<int, D> coord;
				int                idx = 0;
				for (size_t x = 0; x < D; x++) {
					coord[x] = (i / strides[x]) % n;
					idx += (coord[x] + starts[x]) / 2 * strides[x];
				}
				u_fine[fine_idx + i] += u_coarse[coarse_idx + idx];
			}
		} else {
			for (int i = 0; i < npow[D]; i++) {
				u_fine[fine_idx + i] += u_coarse[coarse_idx + i];
			}
		}
	}
	VecRestoreArray(fine, &u_fine);
	VecRestoreArray(coarse_tmp, &u_coarse);
}
} // namespace GMG
#endif
