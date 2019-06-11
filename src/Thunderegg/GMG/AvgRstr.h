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

#ifndef GMGAVGRSTR_H
#define GMGAVGRSTR_H
#include <Thunderegg/DomainCollection.h>
#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/Restrictor.h>
#include <memory>
namespace GMG
{
/**
 * @brief Restrictor that averages the corresponding fine cells into each coarse cell.
 */
template <size_t D> class AvgRstr : public Restrictor<D>
{
	private:
	/**
	 * @brief The coarse DomainCollection that is being restricted to.
	 */
	std::shared_ptr<DomainCollection<D>> coarse_dc;
	/**
	 * @brief The fine DomainCollection that is being restricted from.
	 */
	std::shared_ptr<DomainCollection<D>> fine_dc;
	/**
	 * @brief The communication package for restricting between levels.
	 */
	std::shared_ptr<InterLevelComm<D>> ilc;

	public:
	/**
	 * @brief Create new AvgRstr object.
	 *
	 * @param coarse_dc the DomainColleciton that is being restricted to.
	 * @param fine_dc the DomainCollection that is being restricted from.
	 * @param ilc the communcation package for these two levels.
	 */
	AvgRstr(std::shared_ptr<DomainCollection<D>> coarse_dc,
	        std::shared_ptr<DomainCollection<D>> fine_dc, std::shared_ptr<InterLevelComm<D>> ilc);
	/**
	 * @brief restriction function
	 *
	 * @param coarse the output vector that is restricted to.
	 * @param fine the input vector that is restricted.
	 */
	void restrict(std::shared_ptr<Vector<D>> coarse, std::shared_ptr<const Vector<D>> fine) const;
};

template <size_t D>
inline AvgRstr<D>::AvgRstr(std::shared_ptr<DomainCollection<D>> coarse_dc,
                           std::shared_ptr<DomainCollection<D>> fine_dc,
                           std::shared_ptr<InterLevelComm<D>>   ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}
template <size_t D>
inline void AvgRstr<D>::restrict(std::shared_ptr<Vector<D>>       coarse,
                                 std::shared_ptr<const Vector<D>> fine) const
{
	std::shared_ptr<Vector<D>> coarse_local = ilc->getNewCoarseDistVec();

	for (ILCFineToCoarseMetadata<D> data : ilc->getFineDomains()) {
		PatchInfo<D> &pinfo             = *data.pinfo;
		LocalData<D>  coarse_local_data = coarse_local->getLocalData(data.local_index);
		LocalData<D>  fine_data         = fine->getLocalData(pinfo.local_index);

		if (pinfo.hasCoarseParent()) {
			Orthant<D>         orth = pinfo.orth_on_parent;
			std::array<int, D> starts;
			for (size_t i = 0; i < D; i++) {
				starts[i] = orth.isOnSide(2 * i) ? 0 : coarse_local_data.getLengths()[i];
			}

			nested_loop<D>(fine_data.getStart(), fine_data.getEnd(),
			               [&](const std::array<int, D> &coord) {
				               std::array<int, D> coarse_coord;
				               for (size_t x = 0; x < D; x++) {
					               coarse_coord[x] = (coord[x] + starts[x]) / 2;
				               }
				               coarse_local_data[coarse_coord] += fine_data[coord] / (1 << D);
			               });
		} else {
			nested_loop<D>(
			fine_data.getStart(), fine_data.getEnd(),
			[&](const std::array<int, D> &coord) { coarse_local_data[coord] += fine_data[coord]; });
		}
	}

	// scatter
	coarse->set(0);
	ilc->scatterReverse(coarse_local, coarse);
}
} // namespace GMG
#endif
