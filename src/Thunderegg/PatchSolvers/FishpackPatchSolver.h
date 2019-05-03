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

#ifndef FISHPACKPATCHSOLVER_H
#define FISHPACKPATCHSOLVER_H
#include "PatchSolvers/PatchSolver.h"
class FishpackPatchSolver : public PatchSolver<2>
{
	double lambda = 0;

	public:
	FishpackPatchSolver(double lambda = 0) { this->lambda = lambda; }
	~FishpackPatchSolver() {}
	void addDomain(SchurDomain<2> &d) {}
	void domainSolve(std::deque<SchurDomain<2>> &domains, const Vec f, Vec u, const Vec gamma)
	{
		for (SchurDomain<2> &d : domains) {
			solve(d, f, u, gamma);
		}
	}
	void solve(SchurDomain<2> &d, const Vec f, Vec u, const Vec gamma);
};
#endif
