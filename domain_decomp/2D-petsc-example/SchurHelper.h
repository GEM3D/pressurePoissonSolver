#ifndef SCHURHELPER_H
#define SCHURHELPER_H
#include "DomainCollection.h"
#include "Interpolator.h"
#include "MyTypeDefs.h"
#include "PatchOperator.h"
#include "PatchSolvers/PatchSolver.h"
#include <petscmat.h>
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
	DomainCollection dc;

	std::shared_ptr<Vec> local_gamma;
	std::shared_ptr<Vec> local_interp;
	/**
	 * @brief The MPI communicator used.
	 */
	Teuchos::RCP<const Teuchos::Comm<int>> comm;

	/**
	 * @brief Interpolates to interface values
	 */
	Teuchos::RCP<Interpolator> interpolator;

	/**
	 * @brief The patch operator
	 */
	Teuchos::RCP<PatchOperator> op;

	/**
	 * @brief The patch solver
	 */
	Teuchos::RCP<PatchSolver> solver;

	typedef std::function<void(int, int, Teuchos::RCP<std::valarray<double>>, bool, bool)> inserter;
	void assembleMatrix(inserter insertBlock);

	public:
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param dc the DomainCollection
	 * @param comm the teuchos communicator
	 */
	SchurHelper(DomainCollection dc, Teuchos::RCP<const Teuchos::Comm<int>> comm,
	            Teuchos::RCP<PatchSolver> solver, Teuchos::RCP<PatchOperator> op,
	            Teuchos::RCP<Interpolator> interpolator);

	/**
	 * @brief Solve with a given set of interface values
	 *
	 * @param f the rhs vector
	 * @param u the vector to put solution in
	 * @param gamma the interface values to use
	 * @param diff the resulting difference
	 */
	void solveWithInterface(const Vec f, Vec u, const Vec gamma, Vec diff);

	/**
	 * @brief Apply patch operator with a given set of interface values
	 *
	 * @param u the solution vector to use
	 * @param gamma the interface values to use
	 * @param f the resulting rhs vector
	 */
	void applyWithInterface(const Vec u, const Vec gamma, Vec f);

	/**
	 * @brief Form the Schur complement matrix
	 *
	 * @return the formed matrix
	 */
	std::shared_ptr<Mat> formCRSMatrix();
};
#endif
