#ifndef IFACEMATRIXHELPER_H
#define IFACEMATRIXHELPER_H
#include "DomainCollection.h"
#include <petscmat.h>
class IfaceMatrixHelper
{
	private:
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
	IfaceMatrixHelper(DomainCollection dc);

	/**
	 * @brief Form the matrix
	 *
	 * @return the formed matrix
	 */
	PW_explicit<Mat> formCRSMatrix(double lambda = 0);
};
#endif
