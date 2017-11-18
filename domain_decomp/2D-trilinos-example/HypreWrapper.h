#ifndef HYPREWRAPPER_H
#define HYPREWRAPPER_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

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
	HypreWrapper(Teuchos::RCP<matrix_type> A, const DomainSignatureCollection &dsc, int n,
	             double tol);
	~HypreWrapper();
	void solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b);
};
#endif
