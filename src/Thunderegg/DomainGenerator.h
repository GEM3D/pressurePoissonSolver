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

#ifndef THUNDEREGG_DOMAINGENERATOR_H
#define THUNDEREGG_DOMAINGENERATOR_H
#include <Thunderegg/Domain.h>
/**
 * @brief Generates Domain objects.
 *
 * This class is intended to wrap around octree/quadtree libraries and provide Thunderegg with the
 * necessary patch information. See P4estDG for the implimentation with the p4est library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class DomainGenerator
{
	public:
	/**
	 * @brief Destroy the DomainGenerator object
	 */
	~DomainGenerator(){};
	/**
	 * @brief Return the finest domain
	 */
	virtual std::shared_ptr<Domain<D>> getFinestDomain() = 0;
	/**
	 * @brief return true if there is a coarser domain to be generated.
	 */
	virtual bool hasCoarserDomain() = 0;
	/**
	 * @brief Return a new coarser domain
	 */
	virtual std::shared_ptr<Domain<D>> getCoarserDomain() = 0;
};
#endif
