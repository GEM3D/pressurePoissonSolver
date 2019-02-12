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

#ifndef GMGHELPER_H
#define GMGHELPER_H
#include "Cycle.h"
#include "DomainCollection.h"
#include "SchurHelper.h"
#include <petscpc.h>
namespace GMG
{
class Helper2d
{
	private:
	std::unique_ptr<Cycle> cycle;

	void apply(Vec f, Vec u);

	public:
	static int multiply(PC A, Vec f, Vec u)
	{
		Helper2d *gh = nullptr;
		PCShellGetContext(A, (void **) &gh);
		VecScale(u, 0);
		gh->apply(f, u);
		return 0;
	}

	Helper2d(int n, std::vector<std::shared_ptr<DomainCollection<2>>> domains,
	         std::shared_ptr<SchurHelper<2>> sh, std::string config_file);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
} // namespace GMG
#endif
