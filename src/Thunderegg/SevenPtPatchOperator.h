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

#ifndef SEVENPTPATCHOPERATOR_H
#define SEVENPTPATCHOPERATOR_H
#include <Thunderegg/PatchOperator.h>
class SevenPtPatchOperator : public PatchOperator<3>
{
	public:
	void apply(SchurInfo<3> &d, std::shared_ptr<const Vector<3>> u,
	           std::shared_ptr<const Vector<2>> gamma, std::shared_ptr<Vector<3>> f);
};
#endif
