#ifndef GMGInterpolator_H
#define GMGInterpolator_H
#include "PW.h"
#include "petscvec.h"
namespace GMG
{
/**
 * @brief base class for interpolation operators from finer levels to coarser levels.
 */
class Interpolator
{
	public:
	/**
	 * @brief Virtual interpolation operation that needs to be implemented in derived classes.
	 *
	 * @param coarse the input vector from the coarser level.
	 * @param fine the output vector for the fine level.
	 */
	virtual void interpolate(PW<Vec> coarse, PW<Vec> fine) const = 0;
};
} // namespace GMG
#endif
