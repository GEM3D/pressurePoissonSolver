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

#ifndef THUNDEREGG_P4ESTDOMGEN_H
#define THUNDEREGG_P4ESTDOMGEN_H
#include <Thunderegg/DomainGenerator.h>
#include <functional>
#include <list>
#include <p4est_extended.h>
/**
 * @brief Generates Domain objects form a given p4est object
 */
class P4estDomGen : public DomainGenerator<2>
{
	public:
	/**
	 * @brief Maps coordinate in a block to a coordinate in the domain.
	 *
	 * Each block is treated as a unit square. The input wil be the block number and a coordinate
	 * withing that unit square.
	 *
	 * @param block_no the block number
	 * @param unit_x the x coordinate in the block
	 * @param unit_y the y coordinate in the block
	 * @param x the resulting x coordinate of the mapping function
	 * @param y the resulting y coordinate of the mapping function
	 */
	using BlockMapFunc
	= std::function<void(int block_no, double unit_x, double unit_y, double &x, double &y)>;

	private:
	/**
	 * @brief copy of p4est tree
	 */
	p4est_t *my_p4est;
	/**
	 * @brief List of the domains
	 *
	 * Finest domain is stored in front
	 */
	std::list<std::shared_ptr<Domain<2>>> domain_list;

	/**
	 * @brief The current level that has been generated.
	 *
	 * Will start with num_levels-1
	 */
	int curr_level;
	/**
	 * @brief The number of levels the p4est quad-tree
	 */
	int num_levels;
	/**
	 * @brief The dimensions of each patch
	 */
	std::array<int, 2> ns;
	/**
	 * @brief The length of a block on the x-axis
	 */
	double x_scale;
	/**
	 * @brief The length of a block on the y-axis
	 */
	double y_scale;
	/**
	 * @brief The block Mapping function being used.
	 */
	BlockMapFunc bmf;
	/**
	 * @brief The function used to set neumann boundary conditions
	 */
	IsNeumannFunc<2> inf;

	/**
	 * @brief Get a new coarser level and add it to the end of domain_list
	 */
	void extractLevel();

	public:
	/**
	 * @brief Construct a new P4estDomGen object
	 *
	 * @param p4est the p4est object
	 * @param ns the number of cells in each direction
	 * @param inf the function used to set neumann boundary conditions
	 * @param bmf the function used to map the blocks to the domain
	 */
	P4estDomGen(p4est_t *p4est, const std::array<int, 2> &ns, IsNeumannFunc<2> inf,
	            BlockMapFunc bmf);
	~P4estDomGen();
	std::shared_ptr<Domain<2>> getFinestDomain();
	bool                       hasCoarserDomain();
	std::shared_ptr<Domain<2>> getCoarserDomain();
};
#endif
