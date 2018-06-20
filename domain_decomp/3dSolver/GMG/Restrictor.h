#ifndef GMGRestrictor_H
#define GMGRestrictor_H
#include "PW.h"
#include "petscvec.h"
namespace GMG
{
/**
 * @brief Base class for multi-grid restriction operators.
 */
class Restrictor
{
	public:
	/**
	 * @brief Virtual function that base classes have to implement.
	 *
	 * @param coarse the output vector that is restricted to.
	 * @param fine the input vector that is restricted.
	 */
	virtual void restrict(PW<Vec> coarse, PW<Vec> fine) const = 0;
};
} // namespace GMG
#endif
