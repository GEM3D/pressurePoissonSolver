#ifndef GMGMatOp_H
#define GMGMatOp_H
#include "Operator.h"
#include "petscmat.h"
#include <memory>
namespace GMG
{
/**
 * @brief Wrapper for PETSc Matrix
 */
class MatOp : public Operator
{
	private:
	/**
	 * @brief PETSc Matrix object
	 */
	PW<Mat> matrix;

	public:
	/**
	 * @brief Crate new MatOp
	 *
	 * @param matrix the PETSc matrix
	 */
	MatOp(PW<Mat> matrix)
	{
		this->matrix = matrix;
	}
	/**
	 * @brief Perform matrix/vector multiply.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(PW<Vec> x, PW<Vec> b) const
	{
		MatMult(matrix, x, b);
	}
};
} // namespace GMG
#endif
