#ifndef GMGFFTBlockJacobiSmoother_H
#define GMGFFTBlockJacobiSmoother_H
#include "PW.h"
#include "SchurHelper.h"
#include "petscvec.h"
namespace GMG
{
/**
 * @brief A block Jacobi smoother that uses FFTW solves on each patch. Implemented using the
 * SchurHelper class.
 */
template <size_t D> class FFTBlockJacobiSmoother : public Smoother
{
	private:
	/**
	 * @brief point to the SchurHelper object.
	 */
	std::shared_ptr<SchurHelper<D>> sh;

	public:
	/**
	 * @brief Create new smoother with SchurHelper object
	 *
	 * @param sh pointer to the SchurHelper object
	 */
	FFTBlockJacobiSmoother(std::shared_ptr<SchurHelper<D>> sh)
	{
		this->sh = sh;
	}
	/**
	 * @brief Run an iteration of smoothing.
	 *
	 * @param f the RHS vector
	 * @param u the solution vector, updated upon return.
	 */
	void smooth(PW<Vec> f, PW<Vec> u) const
	{
		sh->solveWithSolution(f, u);
	}
};
} // namespace GMG
#endif
