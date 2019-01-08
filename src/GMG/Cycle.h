#ifndef GMGCycle_H
#define GMGCycle_H
#include "Level.h"
#include <list>
namespace GMG
{
/**
 * @brief Base class for cycles. Includes functions for preparing vectors for finer and coarser
 * levels, and a function to run an iteration of smoothing on a level. Derived cycle classes
 * need to implement the visit function.
 */
class Cycle
{
	private:
	/**
	 * @brief stack of LHS vectors in use. Coarsest at beginning, finest at end.
	 */
	std::list<PW<Vec>> u_vectors;
	/**
	 * @brief stack of RHS vectors in use. Coarsest at beginning, finest at end.
	 */
	std::list<PW<Vec>> f_vectors;
	/**
	 * @brief pointer to the finest level
	 */
	std::shared_ptr<Level> finest_level;

	protected:
	/**
	 * @brief Prepare vectors for coarser level.
	 *
	 * @param level the current level
	 */
	void prepCoarser(const Level &level);
	/**
	 * @brief Prepare vectors for finer level
	 *
	 * @param level the current level
	 */
	void prepFiner(const Level &level);

	/**
	 * @brief run iteration of smoother on solution
	 *
	 * @param level the current level
	 */
	void smooth(const Level &level);

	/**
	 * @brief Virtual visit function that needs to be implemented in derived classes.
	 *
	 * @param level the level currently begin visited.
	 */
	virtual void visit(const Level &level) = 0;

	public:
	/**
	 * @brief Create new cycle object.
	 *
	 * @param finest_level pointer to the finest level object.
	 */
	Cycle(std::shared_ptr<Level> finest_level)
	{
		this->finest_level = finest_level;
	}
	/**
	 * @brief Run one iteration of the cycle.
	 *
	 * @param f the RHS vector.
	 * @param u the current solution vector. Output will be updated solution vector.
	 */
	void apply(Vec f, Vec u)
	{
		f_vectors.push_back(PW<Vec>(f, false));
		u_vectors.push_back(PW<Vec>(u, false));
		visit(*finest_level);
		f_vectors.pop_front();
		u_vectors.pop_front();
	}
};
} // namespace GMG
#endif
