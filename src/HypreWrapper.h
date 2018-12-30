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
