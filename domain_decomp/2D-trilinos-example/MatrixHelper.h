#ifndef MATRIXHELPER_H
#define MATRIXHELPER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
class MatrixHelper
{
	private:
	/**
	 * @brief The MPI communicator used.
	 */
	Teuchos::RCP<const Teuchos::Comm<int>> comm;

	DomainCollection dc;

	public:
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param dc the DomainCollection
	 * @param n number of cells in each direction for each domain
	 * @param h_x the x spacing
	 * @param h_y the y spacing
	 * @param comm the teuchos communicator
	 */
	MatrixHelper(DomainCollection dc, Teuchos::RCP<const Teuchos::Comm<int>> comm);

	/**
	 * @brief Form the matrix
	 *
	 * @return the formed matrix
	 */
	Teuchos::RCP<matrix_type> formCRSMatrix();
};
#endif
