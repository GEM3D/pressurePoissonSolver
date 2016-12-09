#ifndef TRIDIAGSOLVER_H_SEEN
#define TRIDIAGSOLVER_H_SEEN

#include "Domain.h"

class TriDiagSolver
{
	private:
	double left_bnd;
	double right_bnd;

	public:
	/**
	 * @brief create a new solver, with the given boundary conditions
	 *
	 * @param left_boundary the boundary condition on the left end
	 * @param right_boundary the boundary condition on the right end
	 */
	TriDiagSolver(double left_boundary, double right_boundary);

	/**
	 * @brief solve for the given grid
	 *
	 * @param dom the domain to solve over
	 */
	void solve(Domain &grid);
};

#endif
