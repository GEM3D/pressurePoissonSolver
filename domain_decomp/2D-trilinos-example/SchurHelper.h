#ifndef SHURHELPER_H
#define SHURHELPER_H
#include "DomainCollection.h"
#include "Interpolator.h"
#include "MyTypeDefs.h"
#include "PatchOperator.h"
#include "PatchSolvers/PatchSolver.h"
#include <valarray>
/**
 * @brief This class represents a collection of domains that a single processor owns.
 *
 * The purposes of this class:
 *   - Provide a member function for solving with a given interface vector.
 *   - Handle the initialization of the domains.
 *   - Provide member functions for calculating error, residual, etc.
 *   - Provide member functions that generate the Schur complement matrix.
 */
class SchurHelper
{
	private:
	Teuchos::RCP<PatchSolver> solver;
	/**
	 * @brief The MPI communicator used.
	 */
	Teuchos::RCP<const Teuchos::Comm<int>> comm;
	int                                    num_cols = 0;

	typedef std::function<void(int, int, Teuchos::RCP<std::valarray<double>>, bool, bool)> inserter;
	void assembleMatrix(inserter insertBlock);

	public:
	Teuchos::RCP<Interpolator> interpolator;

	Teuchos::RCP<PatchOperator> op;

	void setPatchSolver(Teuchos::RCP<PatchSolver> psolver);

	DomainCollection dc;

	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param dc the DomainCollection
	 * @param n number of cells in each direction for each domain
	 * @param h_x the x spacing
	 * @param h_y the y spacing
	 * @param comm the teuchos communicator
	 */
	SchurHelper(DomainCollection dc, Teuchos::RCP<const Teuchos::Comm<int>> comm);

	/**
	 * @brief Solve with a given set of interface values
	 *
	 * @param gamma the interface values to use
	 * @param diff the resulting difference
	 */
	void solveWithInterface(const vector_type &f, vector_type &u, const vector_type &gamma,
	                        vector_type &diff);

	void applyWithInterface(const vector_type &u, const vector_type &gamma, vector_type &f);

	/**
	 * @brief Form the Schur complement matrix
	 *
	 * @param map the map used in the matrix
	 *
	 * @return the formed matrix
	 */
	void formCRSMatrix(Teuchos::RCP<map_type> map, Teuchos::RCP<matrix_type> &A);
};
#endif
