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

#ifndef GMGVCycle_H
#define GMGVCycle_H
#include <Thunderegg/GMG/Cycle.h>
#include <Thunderegg/tpl/json.hpp>
namespace GMG
{
/**
 * @brief Implementation of a V-cycle
 */
template <size_t D> class VCycle : public Cycle<D>
{
	private:
	int num_pre_sweeps    = 1;
	int num_post_sweeps   = 1;
	int num_coarse_sweeps = 1;

	protected:
	/**
	 * @brief Implements V-cycle. Pre-smooth, visit coarser level and then post-smooth.
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
			for (int i = 0; i < num_post_sweeps; i++) {
				this->smooth(level, u_vectors, f_vectors);
			}
		}
		if (!level.finest()) { this->prepFiner(level, u_vectors, f_vectors); }
	}

	public:
	/**
	 * @brief Create new V-cycle
	 *
	 * @param finest_level a pointer to the finest level
	 */
	VCycle(std::shared_ptr<Level<D>> finest_level, nlohmann::json config_j) : Cycle<D>(finest_level)
	{
		try {
			num_pre_sweeps = config_j.at("pre_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
			throw 0;
		}
		try {
			num_post_sweeps = config_j.at("post_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
			throw 0;
		}
		try {
			num_coarse_sweeps = config_j.at("coarse_sweeps");
		} catch (nlohmann::detail::out_of_range oor) {
			throw 0;
		}
	}
};
} // namespace GMG
#endif
