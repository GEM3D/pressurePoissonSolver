#include "AmgxWrapper.h"
#include <map>
#include <set>
#include <vector>
using namespace std;
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

AmgxWrapper::AmgxWrapper(Teuchos::RCP<matrix_type> A, const DomainSignatureCollection &dsc, int n)
{
	AmgxCrs Acrs;
	Acrs.nnz = A->getNodeNumEntries();
	num_rows = A->getNodeNumRows();
	num_cols = A->getNodeNumCols();
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
     AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));
	MPI_Comm AMGX_MPI_Comm = MPI_COMM_WORLD;
	// app must know how to provide a mapping
	int devices[]   = {/*get_device_id_for_this_rank()*/ 0};
	int num_devices = 1;
	int rank        = 0;
	AMGX_resources_create(&rsrc, cfg, &AMGX_MPI_Comm, 1, &rank);
	AMGX_matrix_create(&gA, rsrc, mode);
	AMGX_vector_create(&gx, rsrc, mode);
	AMGX_vector_create(&gb, rsrc, mode);
	AMGX_solver_create(&solver, rsrc, mode, cfg);

	// create communication pattern
	AmgxMap am(dsc.amgxmap, n);

	AMGX_matrix_comm_from_maps_one_ring(gA, 1, am.num_neighbors, &am.neighbors[0],
	                                          &am.send_sizes[0], &am.send_maps[0], &am.recv_sizes[0],
	                                          &am.recv_maps[0]);
	AMGX_vector_bind(gx, gA);
	AMGX_vector_bind(gb, gA);

	AMGX_matrix_upload_all(gA, num_rows, Acrs.nnz, 1, 1, &Acrs.row_ptrs[0],
	                                      &Acrs.cols[0], (void *) &Acrs.data[0], nullptr);
	AMGX_solver_setup(solver, gA);
}
void AmgxWrapper::solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b)
{
	AMGX_vector_upload(gb, num_rows, 1, (void *) b->get1dViewNonConst().get());
	AMGX_vector_upload(gx, num_rows, 1, (void *) x->get1dViewNonConst().get());
	AMGX_solver_solve_with_0_initial_guess(solver, gb, gx);
	AMGX_vector_download(gx, (void *) x->get1dViewNonConst().get());
}
