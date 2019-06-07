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

#ifndef P4ESTDCG_H
#define P4ESTDCG_H
#include <Thunderegg/DomainCollectionGenerator.h>
#include <functional>
#include <list>
#include <p4est_extended.h>
class P4estDCG : public DomainCollectionGenerator<2>
{
	public:
	using BlockMapFunc = std::function<void(int, double, double, double &, double &)>;

	private:
	bool neumann;
	/**
	 * @brief copy of p4est tree
	 */
	p4est_t *my_p4est;
	/**
	 * @brief finest is stored in front
	 */
	std::list<std::shared_ptr<DomainCollection<2>>> dc_list;

	int                curr_level;
	int                num_levels;
	std::array<int, 2> ns;
	double             x_scale;
	double             y_scale;
	BlockMapFunc       bmf;

	void extractLevel();

	public:
	P4estDCG(p4est_t *p4est, const std::array<int, 2> &ns, bool neumann,
	         BlockMapFunc bmf);
	~P4estDCG();
	std::shared_ptr<DomainCollection<2>> getFinestDC();
	bool                                 hasCoarserDC();
	std::shared_ptr<DomainCollection<2>> getCoarserDC();
};
#endif
