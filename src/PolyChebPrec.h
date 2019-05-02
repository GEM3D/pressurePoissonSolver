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

#ifndef POLYCHEBPREC_H
#define POLYCHEBPREC_H
#include <DomainCollection.h>
#include <MatrixHelper.h>
#include <Operators/Operator.h>
#include <SchurHelper.h>
#include <petscpc.h>
class PolyChebPrec : public Operator<2>
{
	private:
	std::shared_ptr<SchurHelper<3>>      sh;
	std::shared_ptr<DomainCollection<3>> dc;
	double                               interval = 0.95;
	std::vector<double>                  coeffs
	= {4.472135954953655e+00, 5.675247900481234e+00, 3.601012922685066e+00, 2.284885928634731e+00,
	   1.449787551186771e+00, 9.199076055378766e-01, 5.836924189936992e-01, 3.703598469934007e-01,
	   2.349977690621489e-01, 1.491089055767314e-01, 9.461139059090561e-02, 6.003206306517687e-02,
	   3.809106471898141e-02, 2.416923786484517e-02, 1.533567161022980e-02, 1.628851184599676e-02};

	public:
	PolyChebPrec(std::shared_ptr<DomainCollection<3>> dc, std::shared_ptr<SchurHelper<3>> sh);

	void apply(std::shared_ptr<const Vector<2>> x, std::shared_ptr<Vector<2>> b) const;
};
#endif
