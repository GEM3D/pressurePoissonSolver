#ifndef AMGXWRAPPER_H
#define AMGXWRAPPER_H
#include "MyTypeDefs.h"
#include "amgx_c.h"
class AmgxWrapper
{
	private:
	static void print_callback(const char *msg, int length) { std::cout << msg; }
	// library handles
	AMGX_Mode             mode = AMGX_mode_dDDI;
	AMGX_config_handle    cfg;
	AMGX_resources_handle rsrc;
	AMGX_matrix_handle    gA;
	AMGX_vector_handle    gb, gx;
	AMGX_solver_handle    solver;
	// status handling
	AMGX_SOLVE_STATUS status;
    int num_rows;

	public:
	AmgxWrapper(Teuchos::RCP<matrix_type> A);
	void solve(Teuchos::RCP<vector_type> x, Teuchos::RCP<vector_type> b);
};
#endif
