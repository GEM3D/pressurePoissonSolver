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

#ifndef HYPREWRAPPER_H
#define HYPREWRAPPER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <Teuchos_RCP.hpp>
class BoomerPrec : public Tpetra::Operator<scalar_type>
{
	private:
	Teuchos::RCP<matrix_type> A;
	HYPRE_IJMatrix            Aij;
	HYPRE_IJVector            bij;
	HYPRE_IJVector            xij;
	HYPRE_ParCSRMatrix        par_A;
	HYPRE_ParVector           par_b;
	HYPRE_ParVector           par_x;
	HYPRE_Solver              precond;

	int num_rows;
	int num_cols;

	public:
	BoomerPrec(Teuchos::RCP<matrix_type> A, const DomainCollection &dsc, int n, double tol,
	           bool schur);
	~BoomerPrec();
	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	           scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const;
	Teuchos::RCP<const map_type> getDomainMap() const { return A->getDomainMap(); }
	Teuchos::RCP<const map_type> getRangeMap() const { return A->getRangeMap(); }
};
class HypreWrapper
{
	private:
	HYPRE_IJMatrix     Aij;
	HYPRE_IJVector     bij;
	HYPRE_IJVector     xij;
	HYPRE_ParCSRMatrix par_A;
	HYPRE_ParVector    par_b;
	HYPRE_ParVector    par_x;
	HYPRE_Solver       solver;
	HYPRE_Solver       precond;

	int num_rows;
	int num_cols;

	public:
	HypreWrapper(Teuchos::RCP<matrix_type> A, const DomainCollection &dsc, int n, double tol,
	             bool schur);
	~HypreWrapper();
	void solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b);
};
#endif
