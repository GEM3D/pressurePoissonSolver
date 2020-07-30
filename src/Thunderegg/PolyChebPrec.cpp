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

#include "PolyChebPrec.h"
#include <iostream>
using namespace std;
PolyChebPrec::PolyChebPrec(std::shared_ptr<Domain<3>> domain, std::shared_ptr<SchurHelper<3>> sh)
{
	this->sh     = sh;
	this->domain = domain;
}
void PolyChebPrec::apply(std::shared_ptr<const Vector<2>> x, std::shared_ptr<Vector<2>> b) const
{
	std::shared_ptr<PetscVector<3>> f = domain->getNewDomainVec();
	std::shared_ptr<PetscVector<3>> u = domain->getNewDomainVec();

	std::shared_ptr<PetscVector<2>> bk  = sh->getNewSchurVec();
	std::shared_ptr<PetscVector<2>> bk1 = sh->getNewSchurVec();
	std::shared_ptr<PetscVector<2>> bk2 = sh->getNewSchurVec();

	for (int i = coeffs.size() - 1; i > 0; i--) {
		sh->solveAndInterpolateWithInterface(f, u, bk1, bk);
		bk->scaleThenAddScaled(4 / interval, -2, bk1);
		bk->addScaled(coeffs[i], x, -1, bk2);
		auto tmp = bk2;
		bk2      = bk1;
		bk1      = bk;
		bk       = tmp;
	}
	sh->solveAndInterpolateWithInterface(f, u, bk1, b);
	b->scaleThenAddScaled(2 / interval, -1, bk1);
	b->addScaled(coeffs[0], x, -1, bk2);
}
