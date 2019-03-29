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

#ifndef GMGRestrictor_H
#define GMGRestrictor_H
#include <PW.h>
#include <petscvec.h>
namespace GMG
{
/**
 * @brief Base class for multi-grid restriction operators.
 */
class Restrictor
{
	public:
	/**
	 * @brief Virtual function that base classes have to implement.
	 *
	 * @param coarse the output vector that is restricted to.
	 * @param fine the input vector that is restricted.
	 */
	virtual void restrict(PW<Vec> coarse, PW<Vec> fine) const = 0;
};
} // namespace GMG
#endif
