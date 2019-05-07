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

#ifndef GMGWCycle_H
#define GMGWCycle_H
#include <Thunderegg/GMG/Cycle.h>
#include <Thunderegg/tpl/json.hpp>
namespace GMG
{
/**
 * @brief Implementation of a W-cycle
 */
template <size_t D> class WCycle : public Cycle<D>
{
	private:
	int num_pre_sweeps    = 1;
	int num_post_sweeps   = 1;
	int num_coarse_sweeps = 1;
	int num_mid_sweeps    = 1;

	protected:
	/**
	 * @brief Implements W-cycle. Pre-smooth, visit coarser level, smooth, visit coarse level, and
	 * then post-smooth.
	 *
	 * @param level the current level that is being visited.
	 */
	void visit(const Level<D> &level, std::list<std::shared_ptr<Vector<D>>> &u_vectors,
	           std::list<std::shared_ptr<const Vector<D>>> &f_vectors) const
	{
		if (level.coarsest()) {
			for (int i = 0; i < num_coarse_sweeps; i++) {
				this->smooth(level, u_vectors, f_vectors);
			}
		} else {
			for (int i = 0; i < num_pre_sweeps; i++) {
				this->smooth(level, u_vectors, f_vectors);
			}
			this->prepCoarser(level, u_vectors, f_vectors);
			this->visit(level.getCoarser(), u_vectors, f_vectors);
			for (int i = 0; i < num_mid_sweeps; i++) {
				this->smooth(level, u_vectors, f_vectors);
			}
			this->prepCoarser(level, u_vectors, f_vectors);
			this->visit(level.getCoarser(), u_vectors, f_vectors);
			for (int i = 0; i < num_post_sweeps; i++) {
				this->smooth(level, u_vectors, f_vectors);
			}
		}
		if (!level.finest()) { this->prepFiner(level, u_vectors, f_vectors); }
	}

	public:
	/**
	 * @brief Create new W-cycle
	 *
	 * @param finest_level a pointer to the finest level
	 */
	WCycle(std::shared_ptr<Level<D>> finest_level, const CycleOpts &opts) : Cycle<D>(finest_level)
	{
		num_pre_sweeps    = opts.pre_sweeps;
		num_post_sweeps   = opts.post_sweeps;
		num_coarse_sweeps = opts.coarse_sweeps;
		num_mid_sweeps    = opts.mid_sweeps;
	}
};
} // namespace GMG
#endif
