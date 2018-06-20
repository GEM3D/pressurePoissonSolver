#ifndef GMGSmoother_H
#define GMGSmoother_H
#include "PW.h"
#include "petscvec.h"
namespace GMG
{
/**
 * @brief Base class for multi-grid smoothing operators.
 */
class Smoother
{
	public:
	/**
	 * @brief Virtual function that derived classes have to implement.
	 *
	 * @param f the RHS vector
	 * @param u the solution vector, updated upon return.
	 */
	virtual void smooth(PW<Vec> f, PW<Vec> u) const = 0;
};
} // namespace GMG
#endif
