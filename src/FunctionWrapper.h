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

#ifndef FUNCTIONWRAPPER_H
#define FUNCTIONWRAPPER_H
#include "SchurHelper.h"
#include <iostream>
#include <petscpc.h>
template <size_t D> class FuncWrap
{
	public:
	PW<Vec>              u;
	PW<Vec>              f;
	SchurHelper<D> *     sh = nullptr;
	DomainCollection<D> *dc = nullptr;
	FuncWrap()              = default;
	FuncWrap(SchurHelper<D> *sh, DomainCollection<D> *dc)
	{
		f        = dc->getNewDomainVec();
		u        = dc->getNewDomainVec();
		this->sh = sh;
		this->dc = dc;
	}
	static int multiply(Mat A, Vec x, Vec y)
	{
		FuncWrap *w = nullptr;
		MatShellGetContext(A, &w);
		w->sh->solveWithInterface(w->f, w->u, x, y);
		return 0;
	}
	static PW_explicit<Mat> getMatrix(SchurHelper<D> *sh, DomainCollection<D> *dc)
	{
		PW<Mat>   A;
		int       M       = sh->getSchurVecGlobalSize();
		int       m       = sh->getSchurVecLocalSize();
		FuncWrap *wrapper = new FuncWrap(sh, dc);
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, wrapper, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) multiply);
		return A;
	}
};
template <size_t D> 
class FullFuncWrap
{
	public:
	SchurHelper<D> *     sh = nullptr;
	DomainCollection<D> *dc = nullptr;
	FullFuncWrap(SchurHelper<D> *sh, DomainCollection<D> *dc)
	{
		this->sh = sh;
		this->dc = dc;
	}
	static int multiply(Mat A, Vec x, Vec y)
	{
		FullFuncWrap *w = nullptr;
		MatShellGetContext(A, &w);
		w->sh->apply(x, y);
		return 0;
	}
	static PW_explicit<Mat> getMatrix(SchurHelper<D> *sh, DomainCollection<D> *dc)
	{
		PW<Mat>       A;
		int           M       = dc->getGlobalNumCells();
		int           m       = dc->getLocalNumCells();
		FullFuncWrap *wrapper = new FullFuncWrap(sh, dc);
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, wrapper, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) multiply);
		return A;
	}
};
class SchwarzPrec
{
	public:
	SchurHelper<3> *     sh = nullptr;
	DomainCollection<3> *dc = nullptr;
	SchwarzPrec()           = default;
	SchwarzPrec(SchurHelper<3> *sh, DomainCollection<3> *dc)
	{
		this->sh = sh;
		this->dc = dc;
	}
	void apply(Vec f, Vec u)
	{
		sh->solveWithSolution(f, u);
	}
	static int multiply(PC A, Vec f, Vec u)
	{
		SchwarzPrec *w = nullptr;
		PCShellGetContext(A, (void **) &w);
		VecScale(u, 0);
		w->apply(f, u);
		return 0;
	}
	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};

#endif
