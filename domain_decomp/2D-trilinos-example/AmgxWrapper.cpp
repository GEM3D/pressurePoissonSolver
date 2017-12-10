#include "AmgxWrapper.h"
#include <map>
#include <set>
#include <vector>
using namespace std;
struct AmgxCrs {
	int                  n_global;
	int                  n;
	int                  nnz;
	int                  block_dimx = 1;
	int                  block_dimy = 1;
	std::vector<int>     row_ptrs;
	std::vector<int64_t> cols;
	std::vector<double>  data;
	double *             diag_data = nullptr;
};

AmgxWrapper::AmgxWrapper(Teuchos::RCP<matrix_type> A, const DomainCollection &dc, int n)
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
		size_t n;
		n                              = A->getNumEntriesInLocalRow(i);
		int                        row = A->getRowMap()->getGlobalElement(i);
		vector<int>                cols(n);
		vector<double>             data(n);
		Teuchos::ArrayView<int>    inds_view(&cols[0], n);
		Teuchos::ArrayView<double> vals_view(&data[0], n);
		A->getGlobalRowCopy(row, inds_view, vals_view, n);
		for (size_t j = 0; j < n; j++) {
			Acrs.data[curr_pos] = data[j];
			Acrs.cols[curr_pos] = cols[j];
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
	AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "communicator=MPI"));
	// app must know how to provide a mapping
	// int devices[]   = {/*get_device_id_for_this_rank()*/ 0};
	// int num_devices = 1;
	int rank = dc.comm->getRank()%2;
MPI_Comm_dup(MPI_COMM_WORLD,&AMGX_MPI_COMM);
	AMGX_resources_create(&rsrc, cfg, &AMGX_MPI_COMM, 1, &rank);
	AMGX_matrix_create(&gA, rsrc, AMGX_mode_dDDI);
	AMGX_vector_create(&gx, rsrc, AMGX_mode_dDDI);
	AMGX_vector_create(&gb, rsrc, AMGX_mode_dDDI);
	AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, cfg);

	int nrings = 0;
	AMGX_config_get_default_number_of_rings(cfg, &nrings);
	int n_global = A->getGlobalNumRows();
	cerr<< "UPLOADING!!!!!!!!!!! "<<n_global<<" " <<num_rows<<endl;
	AMGX_matrix_upload_all_global(gA, n_global, num_rows, Acrs.nnz, 1, 1, &Acrs.row_ptrs[0],
	                              &Acrs.cols[0], (void *) &Acrs.data[0], nullptr, nrings, nrings,
	                              nullptr);

	cerr<< "UPLOADED!!!!!!!!!!!"<<endl;
	cerr<< "BINDING!!!!!!!!!!! "<<endl;
	AMGX_vector_bind(gx, gA);
	AMGX_vector_bind(gb, gA);
	cerr<< "BOUND!!!!!!!!!!! "<<endl;

	AMGX_solver_setup(solver, gA);
}
AmgxWrapper::~AmgxWrapper()
{
	AMGX_solver_destroy(solver);
	AMGX_matrix_destroy(gA);
	AMGX_vector_destroy(gx);
	AMGX_vector_destroy(gb);
	AMGX_resources_destroy(rsrc);
	AMGX_config_destroy(cfg);
	AMGX_finalize_plugins();
	AMGX_finalize();
}
void AmgxWrapper::solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b)
{
	AMGX_vector_upload(gb, num_rows, 1, (void *) b->get1dViewNonConst().get());
	AMGX_vector_upload(gx, num_rows, 1, (void *) x->get1dViewNonConst().get());
	AMGX_solver_solve(solver, gb, gx);
	AMGX_vector_download(gx, (void *) x->get1dViewNonConst().get());
}
