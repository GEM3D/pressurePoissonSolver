/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef GMGCycle_H
#define GMGCycle_H
#include <GMG/Level.h>
#include <list>
#include <PetscVector.h>
namespace GMG
{
/**
 * @brief Base class for cycles. Includes functions for preparing vectors for finer and coarser
 * levels, and a function to run an iteration of smoothing on a level. Derived cycle classes
 * need to implement the visit function.
 */
template <size_t D> class Cycle
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
	std::shared_ptr<Level<D>> finest_level;

	protected:
	/**
	 * @brief Prepare vectors for coarser level.
	 *
	 * @param level the current level
	 */
	void prepCoarser(const Level<D> &level)
	{
		// calculate residual
		PW<Vec> r = level.getVectorGenerator()->getNewVector();
		level.getOperator().apply(u_vectors.front(), r);
		VecAYPX(r, -1, f_vectors.front());
		// create vectors for coarser levels
		PW<Vec> new_u = level.getCoarser().getVectorGenerator()->getNewVector();
		PW<Vec> new_f = level.getCoarser().getVectorGenerator()->getNewVector();
        std::shared_ptr<Vector<D>> r_vec(new PetscVector<D>(r,16));
        std::shared_ptr<Vector<D>> f_vec(new PetscVector<D>(new_f,16));
		level.getRestrictor().restrict(f_vec, r_vec);
		u_vectors.push_front(new_u);
		f_vectors.push_front(new_f);
	}
	/**
	 * @brief Prepare vectors for finer level
	 *
	 * @param level the current level
	 */
	void prepFiner(const Level<D> &level)
	{
		PW<Vec> old_u = u_vectors.front();
		u_vectors.pop_front();
		f_vectors.pop_front();
        std::shared_ptr<Vector<D>> old_u_vec(new PetscVector<D>(old_u,16));
        std::shared_ptr<Vector<D>> u_vec(new PetscVector<D>(u_vectors.front(),16));
		level.getInterpolator().interpolate(old_u_vec, u_vec);
	}

	/**
	 * @brief run iteration of smoother on solution
	 *
	 * @param level the current level
	 */
	void smooth(const Level<D> &level)
	{
		level.getSmoother().smooth(f_vectors.front(), u_vectors.front());
	}

	/**
	 * @brief Virtual visit function that needs to be implemented in derived classes.
	 *
	 * @param level the level currently begin visited.
	 */
	virtual void visit(const Level<D> &level) = 0;

	public:
	/**
	 * @brief Create new cycle object.
	 *
	 * @param finest_level pointer to the finest level object.
	 */
	Cycle(std::shared_ptr<Level<D>> finest_level)
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
