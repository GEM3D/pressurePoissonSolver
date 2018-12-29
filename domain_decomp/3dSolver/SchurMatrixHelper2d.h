#ifndef SCHURMATRIXHELPER2D_H
#define SCHURMATRIXHELPER2D_H
#include "Interpolator.h"
#include "PatchOperator.h"
#include "PatchSolvers/PatchSolver.h"
#include "SchurHelper.h"
#include <memory>
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
class SchurMatrixHelper2d
{
	private:
	std::shared_ptr<SchurHelper<2>> sh;
	int                             n;

	typedef std::function<void(int, int, std::shared_ptr<std::valarray<double>>, bool, bool)>
	     inserter;
	void assembleMatrix(inserter insertBlock);

	public:
	SchurMatrixHelper2d(std::shared_ptr<SchurHelper<2>> sh)
	{
		this->sh = sh;
		n        = sh->getN();
	}
	/**
	 * @brief Form the Schur complement matrix
	 *
	 * @return the formed matrix
	 */
	PW_explicit<Mat> formCRSMatrix();
};
#endif
