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

#ifndef TRILININTERP_H
#define TRILININTERP_H
#include <Thunderegg/IfaceInterp.h>
class TriLinInterp : public IfaceInterp<3>
{
	public:
	void interpolate(const std::vector<SchurInfo<3>> &patches, std::shared_ptr<const Vector<3>> u,
	                 std::shared_ptr<Vector<2>> interp);
	void interpolate(SchurInfo<3> &d, std::shared_ptr<const Vector<3>> u,
	                 std::shared_ptr<Vector<2>> interp);
	void interpolate(SchurInfo<3> &d, Side<3> s, int local_index, IfaceType<3> itype,
	                 std::shared_ptr<const Vector<3>> u, std::shared_ptr<Vector<2>> interp);
};
#endif
