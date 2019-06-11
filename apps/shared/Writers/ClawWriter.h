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

#ifndef CLAWWRITER_H
#define CLAWWRITER_H
#include <Thunderegg/DomainCollection.h>
class ClawWriter
{
	private:
	DomainCollection<2> dc;
	void writePatch(PatchInfo<2> &d, std::ostream &os, double *u_view, double *resid_view);

	public:
	ClawWriter(DomainCollection<2> &dc);
	void write(Vec u, Vec resid);
};
#endif
