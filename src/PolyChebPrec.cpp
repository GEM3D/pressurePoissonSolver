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
PolyChebPrec::PolyChebPrec(SchurHelper<3> &sh, DomainCollection<3> &dc)
{
	this->sh = &sh;
	this->dc = &dc;
}
void PolyChebPrec::apply(Vec b, Vec y)
{
	std::shared_ptr<PetscVector<2>> y_vec(new PetscVector<2>(y, sh->getLengths(), false));

	std::shared_ptr<PetscVector<3>> f = dc->getNewDomainVec();
	std::shared_ptr<PetscVector<3>> u = dc->getNewDomainVec();

	std::shared_ptr<PetscVector<2>> bk  = sh->getNewSchurVec();
	std::shared_ptr<PetscVector<2>> bk1 = sh->getNewSchurVec();
	std::shared_ptr<PetscVector<2>> bk2 = sh->getNewSchurVec();

	for (int i = coeffs.size() - 1; i > 0; i--) {
		sh->solveAndInterpolateWithInterface(f, u, bk1, bk);
		VecAXPBY(bk->vec, -2, 4 / interval, bk1->vec);
		VecAXPBYPCZ(bk->vec, coeffs[i], -1.0, 1.0, b, bk2->vec);
		auto tmp = bk2;
		bk2      = bk1;
		bk1      = bk;
		bk       = tmp;
	}
	sh->solveAndInterpolateWithInterface(f, u, bk1, y_vec);
	VecAXPBY(y, -1, 2 / interval, bk1->vec);
	VecAXPBYPCZ(y, coeffs[0], -1.0, 1.0, b, bk2->vec);
}
