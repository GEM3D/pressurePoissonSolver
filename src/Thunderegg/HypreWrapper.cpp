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

#include "HypreWrapper.h"
#include <HYPRE_parcsr_mv.h>
#include <map>
#include <set>
#include <vector>
using namespace std;
BoomerPrec::BoomerPrec(Teuchos::RCP<matrix_type> A, const DomainCollection &dsc, int n, double tol,
                       bool schur)
{
	this->A    = A;
	int ilower = A->getRowMap()->getMinGlobalIndex();
	int iupper = A->getRowMap()->getMaxGlobalIndex();
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &Aij);
	HYPRE_IJMatrixSetObjectType(Aij, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(Aij);

	num_rows = A->getNodeNumRows();
	for (int i = 0; i < num_rows; i++) {
		size_t n   = A->getNumEntriesInLocalRow(i);
		int    row = A->getRowMap()->getGlobalElement(i);

		vector<int>                cols(n);
		vector<double>             data(n);
		Teuchos::ArrayView<int>    inds_view(&cols[0], n);
		Teuchos::ArrayView<double> vals_view(&data[0], n);
		A->getGlobalRowCopy(row, inds_view, vals_view, n);

		int ncols = n;
		HYPRE_IJMatrixSetValues(Aij, 1, &ncols, &row, &cols[0], &data[0]);
	}
	HYPRE_IJMatrixAssemble(Aij);
	HYPRE_IJMatrixGetObject(Aij, (void **) &par_A);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xij);
	HYPRE_IJVectorSetObjectType(xij, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xij);
	HYPRE_IJVectorAssemble(xij);
	HYPRE_IJVectorGetObject(xij, (void **) &par_x);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bij);
	HYPRE_IJVectorSetObjectType(bij, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bij);
	HYPRE_IJVectorAssemble(bij);
	HYPRE_IJVectorGetObject(bij, (void **) &par_b);

	if (schur) {
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetStrongThreshold(precond, .9);
		HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);
		HYPRE_BoomerAMGSetRelaxType(precond, 0);
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 0, .25);
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetTol(precond, 0.0);
		HYPRE_BoomerAMGSetPrintLevel(precond, 1);
		HYPRE_BoomerAMGSetMaxIter(precond, 1);
	} else {
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
		HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);
		HYPRE_BoomerAMGSetRelaxType(precond, 3);
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 3, 3);
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetTol(precond, 0.0);
		HYPRE_BoomerAMGSetPrintLevel(precond, 1);
		HYPRE_BoomerAMGSetMaxIter(precond, 1);
	}

	// set preconditioner
	HYPRE_BoomerAMGSetup(precond, par_A, par_b, par_x);
}

void BoomerPrec::apply(const vector_type &b, vector_type &x, Teuchos::ETransp mode,
                       scalar_type alpha, scalar_type beta) const
{
	double *xvals = x.get1dViewNonConst().get();
	auto    xinds = x.getMap()->getMyGlobalIndices();
	HYPRE_IJVectorSetValues(xij, num_rows, &xinds[0], xvals);
	HYPRE_IJVectorAssemble(xij);

	const double *bvals = b.get1dView().get();
	auto          binds = b.getMap()->getMyGlobalIndices();
	HYPRE_IJVectorSetValues(bij, num_rows, &binds[0], bvals);
	HYPRE_IJVectorAssemble(bij);

	HYPRE_BoomerAMGSolve(precond, par_A, par_b, par_x);

	HYPRE_IJVectorGetValues(xij, num_rows, &xinds[0], xvals);
}
HypreWrapper::HypreWrapper(Teuchos::RCP<matrix_type> A, const DomainCollection &dsc, int n,
                           double tol, bool schur)
{
	int ilower = A->getRowMap()->getMinGlobalIndex();
	int iupper = A->getRowMap()->getMaxGlobalIndex();
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &Aij);
	HYPRE_IJMatrixSetObjectType(Aij, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(Aij);

	num_rows = A->getNodeNumRows();
	for (int i = 0; i < num_rows; i++) {
		size_t n   = A->getNumEntriesInLocalRow(i);
		int    row = A->getRowMap()->getGlobalElement(i);

		vector<int>                cols(n);
		vector<double>             data(n);
		Teuchos::ArrayView<int>    inds_view(&cols[0], n);
		Teuchos::ArrayView<double> vals_view(&data[0], n);
		A->getGlobalRowCopy(row, inds_view, vals_view, n);

		int ncols = n;
		HYPRE_IJMatrixSetValues(Aij, 1, &ncols, &row, &cols[0], &data[0]);
	}
	HYPRE_IJMatrixAssemble(Aij);
	HYPRE_IJMatrixGetObject(Aij, (void **) &par_A);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xij);
	HYPRE_IJVectorSetObjectType(xij, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xij);
	HYPRE_IJVectorAssemble(xij);
	HYPRE_IJVectorGetObject(xij, (void **) &par_x);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bij);
	HYPRE_IJVectorSetObjectType(bij, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bij);
	HYPRE_IJVectorAssemble(bij);
	HYPRE_IJVectorGetObject(bij, (void **) &par_b);

	if (schur) {
		HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
		HYPRE_ParCSRGMRESSetMaxIter(solver, 1000);
		HYPRE_ParCSRGMRESSetTol(solver, tol);
		HYPRE_ParCSRGMRESSetPrintLevel(solver, 3);
		HYPRE_ParCSRGMRESSetLogging(solver, 1);

		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
		HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);
		HYPRE_BoomerAMGSetRelaxType(precond, 0);
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 0, 3);
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetTol(precond, 0.0);
		HYPRE_BoomerAMGSetPrintLevel(precond, 1);
		HYPRE_BoomerAMGSetMaxIter(precond, 1);
	} else {
		HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
		HYPRE_ParCSRGMRESSetMaxIter(solver, 1000);
		HYPRE_ParCSRGMRESSetTol(solver, tol);
		HYPRE_ParCSRGMRESSetPrintLevel(solver, 3);
		HYPRE_ParCSRGMRESSetLogging(solver, 1);

		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
		HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);
		HYPRE_BoomerAMGSetRelaxType(precond, 3);
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 3, 3);
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetTol(precond, 0.0);
		HYPRE_BoomerAMGSetPrintLevel(precond, 1);
		HYPRE_BoomerAMGSetMaxIter(precond, 1);
	}

	// set preconditioner
	HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);

	HYPRE_ParCSRGMRESSetup(solver, par_A, par_b, par_x);
}
HypreWrapper::~HypreWrapper()
{
	HYPRE_IJMatrixDestroy(Aij);
	HYPRE_IJVectorDestroy(xij);
	HYPRE_IJVectorDestroy(bij);
}
BoomerPrec::~BoomerPrec()
{
	HYPRE_IJMatrixDestroy(Aij);
	HYPRE_IJVectorDestroy(xij);
	HYPRE_IJVectorDestroy(bij);
}
void HypreWrapper::solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b)
{
	double *xvals = x->get1dViewNonConst().get();
	auto    xinds = x->getMap()->getMyGlobalIndices();
	HYPRE_IJVectorSetValues(xij, num_rows, &xinds[0], xvals);
	HYPRE_IJVectorAssemble(xij);

	double *bvals = b->get1dViewNonConst().get();
	auto    binds = b->getMap()->getMyGlobalIndices();
	HYPRE_IJVectorSetValues(bij, num_rows, &binds[0], bvals);
	HYPRE_IJVectorAssemble(bij);

	HYPRE_ParCSRGMRESSolve(solver, par_A, par_b, par_x);

	HYPRE_IJVectorGetValues(xij, num_rows, &xinds[0], xvals);
}
