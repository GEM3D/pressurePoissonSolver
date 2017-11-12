#include "AmgxWrapper.h"
struct AmgxCrs {
	int                 n_global;
	int                 n;
	int                 nnz;
	int                 block_dimx = 1;
	int                 block_dimy = 1;
	std::vector<int>    row_ptrs;
	std::vector<int>    cols;
	std::vector<double> data;
	double *            diag_data = nullptr;
};

AmgxWrapper::AmgxWrapper(Teuchos::RCP<matrix_type> A)
{
	AmgxCrs Acrs;
	Acrs.nnz = A->getNodeNumEntries();
	num_rows = A->getNodeNumRows();
	Acrs.row_ptrs.resize(num_rows + 1);
	Acrs.data.resize(Acrs.nnz);
	Acrs.cols.resize(Acrs.nnz);
	int curr_pos     = 0;
	Acrs.row_ptrs[0] = curr_pos;
	for (int i = 0; i < num_rows; i++) {
		int           n;
		const int *   inds;
		const double *vals;
		A->getLocalRowViewRaw(i, n, inds, vals);
		for (int j = 0; j < n; j++) {
			Acrs.data[curr_pos] = vals[j];
			Acrs.cols[curr_pos] = inds[j];
			curr_pos++;
		}

		Acrs.row_ptrs[i + 1] = curr_pos;
	}
	/* init */
	AMGX_SAFE_CALL(AMGX_initialize());
	AMGX_SAFE_CALL(AMGX_initialize_plugins());
	/* system */
	AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
	AMGX_SAFE_CALL(AMGX_install_signal_handler());
	/* create resources, matrix, vector and solver */
	AMGX_config_create_from_file(&cfg, "amgx.json");
	AMGX_resources_create_simple(&rsrc, cfg);
	AMGX_matrix_create(&gA, rsrc, mode);
	AMGX_vector_create(&gx, rsrc, mode);
	AMGX_vector_create(&gb, rsrc, mode);
	AMGX_solver_create(&solver, rsrc, mode, cfg);

	AMGX_SAFE_CALL(AMGX_matrix_upload_all(gA, num_rows, Acrs.nnz, 1, 1, &Acrs.row_ptrs[0],
	                                      &Acrs.cols[0], (void *) &Acrs.data[0], nullptr));
	AMGX_SAFE_CALL(AMGX_solver_setup(solver, gA));
}
void AmgxWrapper::solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b)
{
	AMGX_SAFE_CALL(AMGX_vector_upload(gb, num_rows, 1, (void *) b->get1dViewNonConst().get()));
	AMGX_SAFE_CALL(AMGX_vector_set_zero(gx, num_rows, 1));
	AMGX_SAFE_CALL(AMGX_solver_solve_with_0_initial_guess(solver, gb, gx));
	AMGX_SAFE_CALL(AMGX_vector_download(gx, (void *) x->get1dViewNonConst().get()));
}
