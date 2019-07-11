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

#ifndef PATCHOPERATOR_H
#define PATCHOPERATOR_H
#include <Thunderegg/SchurInfo.h>
#include <Thunderegg/Vector.h>
template <size_t D> class PatchOperator
{
	public:
	virtual ~PatchOperator() {}
	virtual void applyWithInterface(SchurInfo<D> &d, const LocalData<D> u,
	                                std::shared_ptr<const Vector<D - 1>> gamma, LocalData<D> f)
	= 0;
	virtual void addInterfaceToRHS(SchurInfo<D> &d, std::shared_ptr<const Vector<D - 1>> gamma,
	                               LocalData<D> f)
	= 0;
	virtual void apply(const SchurInfo<D> &d, const LocalData<D> u, LocalData<D> f) = 0;
};
#endif
