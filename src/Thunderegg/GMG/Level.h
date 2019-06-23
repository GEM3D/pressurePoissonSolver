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

#ifndef GMGLEVEL_H
#define GMGLEVEL_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/GMG/Interpolator.h>
#include <Thunderegg/GMG/Restrictor.h>
#include <Thunderegg/GMG/Smoother.h>
#include <Thunderegg/Operators/Operator.h>
#include <Thunderegg/Vector.h>
#include <memory>
namespace GMG
{

/**
 * @brief Represents a level in geometric multi-grid.
 */
template <size_t D> class Level
{
	private:
	/**
	 * @brief the VectorGenerator for this level.
	 */
	std::shared_ptr<VectorGenerator<D>> vg;
	/**
	 * @brief The operator (matrix) for this level.
	 */
	std::shared_ptr<Operator<D>> op;
	/**
	 * @brief The restrictor from this level to the coarser level.
	 */
	std::shared_ptr<Restrictor<D>> restrictor;
	/**
	 * @brief The interpolator from this level to the finer level.
	 */
	std::shared_ptr<Interpolator<D>> interpolator;
	/**
	 * @brief The smoother for this level.
	 */
	std::shared_ptr<Smoother<D>> smoother;
	/**
	 * @brief Pointer to coarser level
	 */
	std::shared_ptr<Level> coarser;
	/**
	 * @brief Pointer to finer level
	 */
	std::weak_ptr<Level> finer;

	public:
	/**
	 * @brief Create a Level object.
	 *
	 * @param dc pointer to the DomainCollection for this level
	 */
	Level(std::shared_ptr<VectorGenerator<D>> vg)
	{
		this->vg = vg;
	}
	/**
	 * @brief Set the restriction operator for restricting from this level to the coarser level.
	 *
	 * @param restrictor the restriction operator.
	 */
	void setRestrictor(std::shared_ptr<Restrictor<D>> restrictor)
	{
		this->restrictor = restrictor;
	}
	/**
	 * @brief Get the restriction operator for this level.
	 *
	 * @return reference to the restrictor
	 */
	const Restrictor<D> &getRestrictor() const
	{
		return *restrictor;
	}
	/**
	 * @brief Set the interpolation operator for interpolating from this level to the finer level.
	 *
	 * @param interpolator the interpolation operator.
	 */
	void setInterpolator(std::shared_ptr<Interpolator<D>> interpolator)
	{
		this->interpolator = interpolator;
	}
	/**
	 * @brief Get the interpolation operator for this level.
	 *
	 * @return Reference to the interpolator.
	 */
	const Interpolator<D> &getInterpolator() const
	{
		return *interpolator;
	}
	/**
	 * @brief Set the operator (matrix) for this level.
	 *
	 * @param op the operator
	 */
	void setOperator(std::shared_ptr<Operator<D>> op)
	{
		this->op = op;
	}
	/**
	 * @brief Get the operator for this level.
	 *
	 * @return Reference to the operator.
	 */
	const Operator<D> &getOperator() const
	{
		return *op;
	}
	/**
	 * @brief Set the smoother for this level.
	 *
	 * @param smoother the smoother
	 */
	void setSmoother(std::shared_ptr<Smoother<D>> smoother)
	{
		this->smoother = smoother;
	}
	/**
	 * @brief Get smoother operator for this level.
	 *
	 * @return Reference to the smoother operator.
	 */
	const Smoother<D> &getSmoother() const
	{
		return *smoother;
	}
	/**
	 * @brief Get DomainCollection for this level.
	 *
	 * @return DomainCollection for this level.
	 */
	const std::shared_ptr<VectorGenerator<D>> &getVectorGenerator() const
	{
		return vg;
	}
	/**
	 * @brief Set pointer to the coarser level.
	 *
	 * @param coarser the pointer to the coarser level.
	 */
	void setCoarser(std::shared_ptr<Level> coarser)
	{
		this->coarser = coarser;
	}
	/**
	 * @brief get reference to the coarser level.
	 *
	 * @return reference to the coarser level.
	 */
	const Level &getCoarser() const
	{
		return *coarser;
	}
	/**
	 * @brief Set the pointer to the finer level.
	 *
	 * @param finer pointer to the finer level.
	 */
	void setFiner(std::shared_ptr<Level> finer)
	{
		this->finer = finer;
	}
	/**
	 * @brief Check if this level is the finest level.
	 *
	 * @return whether or not this level is the finest level.
	 */
	bool finest() const
	{
		return finer.expired();
	}
	/**
	 * @brief Check if this level is the coarsest level.
	 *
	 * @return whether or not this level is the coarsest level.
	 */
	bool coarsest() const
	{
		return coarser == nullptr;
	}
};
} // namespace GMG
#endif
