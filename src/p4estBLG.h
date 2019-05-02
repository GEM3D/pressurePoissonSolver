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

#ifndef P4ESTBLG_H
#define P4ESTBLG_H
#include <Domain.h>
#include <deque>
#include <p4est.h>
#include <p4est_mesh.h>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <set>
#include <vector>
#include <zoltan.h>
class p4estBLG
{
	private:
	void extractLevel(p4est_t *p4est,int level, std::array<int,2> n);

	public:
	using DomainMap = std::map<int, std::shared_ptr<Domain<2>>>;
	std::vector<DomainMap> levels;
	p4estBLG(p4est_t *p4est, std::array<int,2> ns);
};
#endif
