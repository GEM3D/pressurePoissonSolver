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

#include "Cycle.h"
using namespace GMG;
void Cycle::prepCoarser(const Level &level)
{
	// calculate residual
	PW<Vec> r = level.getVectorGenerator()->getNewVector();
	level.getOperator().apply(u_vectors.front(), r);
	VecAYPX(r, -1, f_vectors.front());
	// create vectors for coarser levels
	PW<Vec> new_u = level.getCoarser().getVectorGenerator()->getNewVector();
	PW<Vec> new_f = level.getCoarser().getVectorGenerator()->getNewVector();
	level.getRestrictor().restrict(new_f, r);
	u_vectors.push_front(new_u);
	f_vectors.push_front(new_f);
}
void Cycle::prepFiner(const Level &level)
{
	PW<Vec> old_u = u_vectors.front();
	u_vectors.pop_front();
	f_vectors.pop_front();
	level.getInterpolator().interpolate(old_u, u_vectors.front());
}
void Cycle::smooth(const Level &level)
{
	level.getSmoother().smooth(f_vectors.front(), u_vectors.front());
}
