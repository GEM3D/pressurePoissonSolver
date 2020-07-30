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
#include <Thunderegg/GMG/Level.h>
#include <Thunderegg/Vector.h>
#include <list>
namespace GMG
{
/**
 * @brief Base class for cycles. Includes functions for preparing vectors for finer and coarser
 * levels, and a function to run an iteration of smoothing on a level. Derived cycle classes
 * need to implement the visit function.
 */
template <size_t D> class Cycle : public Operator<D>
{
	private:
	/**
	 * @brief pointer to the finest level
	 */
	std::shared_ptr<Level<D>> finest_level;

	protected:
	using VecList      = std::list<std::shared_ptr<Vector<D>>>;
	using ConstVecList = std::list<std::shared_ptr<const Vector<D>>>;
	/**
	 * @brief stack of LHS vectors in use. Coarsest at beginning, finest at end.
	 */
	/**
	 * @brief stack of RHS vectors in use. Coarsest at beginning, finest at end.
	 */
	/**
	 * @brief Prepare vectors for coarser level.
	 *
	 * @param level the current level
	 */
	void prepCoarser(const Level<D> &level, VecList &u_vectors, ConstVecList &f_vectors) const
	{
		// calculate residual
		std::shared_ptr<Vector<D>> r = level.getVectorGenerator()->getNewVector();
		level.getOperator().apply(u_vectors.front(), r);
		r->scaleThenAdd(-1, f_vectors.front());
		// create vectors for coarser levels
		std::shared_ptr<Vector<D>> new_u = level.getCoarser().getVectorGenerator()->getNewVector();
		std::shared_ptr<Vector<D>> new_f = level.getCoarser().getVectorGenerator()->getNewVector();
		level.getRestrictor().restrict(new_f, r);
		u_vectors.push_front(new_u);
		f_vectors.push_front(new_f);
	}
	/**
	 * @brief Prepare vectors for finer level
	 *
	 * @param level the current level
	 */
	void prepFiner(const Level<D> &level, VecList &u_vectors, ConstVecList &f_vectors) const
	{
		std::shared_ptr<Vector<D>> old_u = u_vectors.front();
		u_vectors.pop_front();
		f_vectors.pop_front();
		level.getInterpolator().interpolate(old_u, u_vectors.front());
	}

	/**
	 * @brief run iteration of smoother on solution
	 *
	 * @param level the current level
	 */
	void smooth(const Level<D> &level, VecList &u_vectors, ConstVecList &f_vectors) const
	{
		level.getSmoother().smooth(f_vectors.front(), u_vectors.front());
	}

	/**
	 * @brief Virtual visit function that needs to be implemented in derived classes.
	 *
	 * @param level the level currently begin visited.
	 */
	virtual void visit(const Level<D> &level, VecList &u_vectors,
	                   ConstVecList &f_vectors) const = 0;

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
	void apply(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u) const
	{
		u->set(0);
		VecList      u_vectors;
		ConstVecList f_vectors;
		f_vectors.push_back(f);
		u_vectors.push_back(u);
		visit(*finest_level, u_vectors, f_vectors);
		f_vectors.pop_front();
		u_vectors.pop_front();
	}
};
} // namespace GMG
#endif
