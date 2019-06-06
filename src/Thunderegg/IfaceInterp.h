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

#ifndef THUNDEREGG_IFACEINTERP_H
#define THUNDEREGG_IFACEINTERP_H
#include <Thunderegg/Iface.h>
#include <Thunderegg/SchurDomain.h>
#include <Thunderegg/Vector.h>
/**
 * @brief An abstract class that interpolates to the interfaces in the Schur compliment system.
 *
 * @tparam D the number of Cartesian dimensions in the patches.
 */
template <size_t D> class IfaceInterp
{
	public:
	virtual ~IfaceInterp() {}

	/**
	 * @brief Given a set of patches, interpolate to their interfaces.
	 *
	 * @param patches the set of patches
	 * @param u the domain vector
	 * @param interp the interface vector
	 */
	virtual void interpolate(const std::vector<SchurDomain<D>> &patches,
	                         std::shared_ptr<const Vector<D>>   u,
	                         std::shared_ptr<Vector<D - 1>>     interp)
	= 0;
	/**
	 * \deprecated
	 * @brief Given a domain vector, interpolate to the interface vector.
	 *
	 * @param d the domain
	 * @param u the domain vector
	 * @param interp the interface vector
	 */
	virtual void interpolate(SchurDomain<D> &d, std::shared_ptr<const Vector<D>> u,
	                         std::shared_ptr<Vector<D - 1>> interp)
	= 0;
	/**
	 * \deprecated
	 * @brief Given a domain vector, interpolate to the interface vector.
	 * @param d
	 * @param s
	 * @param local_index
	 * @param itype
	 * @param u
	 * @param interp
	 */
	virtual void interpolate(SchurDomain<D> &d, Side<D> s, int local_index, IfaceType itype,
	                         std::shared_ptr<const Vector<D>> u,
	                         std::shared_ptr<Vector<D - 1>>   interp)
	= 0;
};
#endif
