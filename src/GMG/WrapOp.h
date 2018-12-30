#ifndef GMGWrapOp_H
#define GMGWrapOp_H
#include "Operator.h"
#include <memory>
namespace GMG
{
/**
 * @brief Wrapper for Matrix free operation
 */
template <size_t D> class WrapOp : public Operator
{
	private:
	/**
	 * @brief PETSc Matrix object
	 */
	std::shared_ptr<SchurHelper<D>> helper;

	public:
	/**
	 * @brief Crate new WrapOp
	 *
	 * @param matrix the PETSc matrix
	 */
	WrapOp(std::shared_ptr<SchurHelper<D>> helper)
	{
		this->helper = helper;
	}
	/**
	 * @brief Perform matrix/vector multiply.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(PW<Vec> x, PW<Vec> b) const
	{
		helper->apply(x, b);
	}
};
} // namespace GMG
#endif
