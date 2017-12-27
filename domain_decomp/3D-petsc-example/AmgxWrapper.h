#ifndef AMGXWRAPPER_H
#define AMGXWRAPPER_H
#include "DomainCollection.h"
#include "amgx_c.h"
#include <petscmat.h>
#include <iostream>
#include <string>
class AmgxWrapper
{
	private:
	static void print_callback(const char *msg, int length)
	{
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0) std::cout << msg;
	}
	// library handles
	MPI_Comm              AMGX_MPI_COMM;
	AMGX_config_handle    cfg;
	AMGX_resources_handle rsrc;
	AMGX_matrix_handle    gA;
	AMGX_vector_handle    gb;
	AMGX_vector_handle    gx;
	AMGX_solver_handle    solver;
	int               num_rows;
        int nrings = 0;

	public:
	AmgxWrapper(std::string filename);
	void setMatrix(Mat A);
	~AmgxWrapper();
	void solve(Vec x, Vec b);
};
#endif
