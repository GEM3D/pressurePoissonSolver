#include "AmgxWrapper.h"
#include <map>
#include <set>
#include <vector>
#include <cuda_runtime.h>
#ifndef NUM_GPU
#define NUM_GPU 1
#endif
using namespace std;
AmgxWrapper::AmgxWrapper(string filename) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());
    /* create resources, matrix, vector and solver */
    AMGX_config_create_from_file(&cfg, filename.c_str());
    AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));
    AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "communicator=MPI"));
    int gpu_id = rank % NUM_GPU;
    AMGX_MPI_COMM=MPI_COMM_WORLD;
    AMGX_resources_create(&rsrc, cfg, &AMGX_MPI_COMM, 1, &gpu_id);
    AMGX_matrix_create(&gA, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&gx, rsrc, AMGX_mode_dDDI);
    AMGX_vector_create(&gb, rsrc, AMGX_mode_dDDI);
    AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, cfg);
    AMGX_config_get_default_number_of_rings(cfg, &nrings);
}
void AmgxWrapper::setMatrix(Mat A) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    PW<Mat> localA;
    MatMPIAIJGetLocalMat(A, MAT_INITIAL_MATRIX, &localA);
    const int *ia;
    const int *ja;
    PetscBool done;
    int nr;
    MatGetRowIJ(localA, 0, PETSC_FALSE, PETSC_FALSE, &nr, &ia, &ja, &done);
    num_rows = nr;
    int nnz = ia[num_rows];


    int m, n;
    MatGetSize(A, &m, &n);

    const int *ranges;
    MatGetOwnershipRanges(A, &ranges);
    vector<int> procs;
    for (int p = 0; p < size; p++) {
        for (int i = ranges[p]; i < ranges[p + 1]; i++) {
            procs.push_back(p);
        }
    }

    vector<double> data(nnz);
    vector<int64_t> cols(nnz);
    vector<int> rows(num_rows + 1);
    int q = 0;
    for (int i = 0; i < num_rows; i++) {
        rows[i] = q;
        int ncols;
        const int *c;
        const double *vals;
        MatGetRow(localA, i, &ncols, &c, &vals);
        for (int j = 0; j < ncols; j++) {
            data[q] = vals[j];
            cols[q] = c[j];
            q++;
        }
    }
    rows[num_rows] = q;


    AMGX_matrix_upload_all_global(gA, n, num_rows, nnz, 1, 1, &rows[0],
                                  &cols[0], (void *)&data[0], nullptr, nrings,
                                  nrings, &procs[0]);

    MatRestoreRowIJ(localA, 0, PETSC_FALSE, PETSC_FALSE, &nr, &ia, &ja, &done);
    AMGX_vector_bind(gx, gA);
    AMGX_vector_bind(gb, gA);

    AMGX_solver_setup(solver, gA);
}
AmgxWrapper::~AmgxWrapper() {
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(gx);
    AMGX_vector_destroy(gb);
    AMGX_matrix_destroy(gA);
    AMGX_resources_destroy(rsrc);
    AMGX_config_destroy(cfg);
    AMGX_reset_signal_handler();
    AMGX_finalize_plugins();
    AMGX_finalize();
}
void AmgxWrapper::solve(Vec x, Vec b) {
    double *x_ptr, *b_ptr;
    VecGetArray(x, &x_ptr);
    VecGetArray(b, &b_ptr);
    AMGX_vector_upload(gb, num_rows, 1, (void *)b_ptr);
    AMGX_vector_upload(gx, num_rows, 1, (void *)x_ptr);
    AMGX_solver_solve_with_0_initial_guess(solver, gb, gx);
    AMGX_vector_download(gx, (void *)x_ptr);
    VecRestoreArray(x, &x_ptr);
    VecRestoreArray(b, &b_ptr);
}
