#ifndef GMGOperator_H
#define GMGOperator_H
#include "PW.h"
#include "petscvec.h"
namespace GMG
{
/**
 * @brief Base class for operators on each level.
 */
class Operator
{
	public:
	/**
	 * @brief Virtual function that base classes have to implement.
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	virtual void apply(PW<Vec> x, PW<Vec> b) const = 0;
};
} // namespace GMG
#endif
