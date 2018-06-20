#ifndef GMGLEVEL_H
#define GMGLEVEL_H
#include "DomainCollection.h"
#include "Interpolator.h"
#include "Operator.h"
#include "Restrictor.h"
#include "Smoother.h"
#include <memory>
namespace GMG
{
/**
 * @brief Represents a level in geometric multi-grid.
 */
class Level
{
	private:
	/**
	 * @brief the DomainCollection for this level.
	 */
	std::shared_ptr<DomainCollection> dc;
	/**
	 * @brief The operator (matrix) for this level.
	 */
	std::shared_ptr<Operator> op;
	/**
	 * @brief The restrictor from this level to the coarser level.
	 */
	std::shared_ptr<Restrictor> restrictor;
	/**
	 * @brief The interpolator from this level to the finer level.
	 */
	std::shared_ptr<Interpolator> interpolator;
	/**
	 * @brief The smoother for this level.
	 */
	std::shared_ptr<Smoother> smoother;
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
	Level(std::shared_ptr<DomainCollection> dc)
	{
		this->dc = dc;
	}
	/**
	 * @brief Set the restriction operator for restricting from this level to the coarser level.
	 *
	 * @param restrictor the restriction operator.
	 */
	void setRestrictor(std::shared_ptr<Restrictor> restrictor)
	{
		this->restrictor = restrictor;
	}
	/**
	 * @brief Get the restriction operator for this level.
	 *
	 * @return reference to the restrictor
	 */
	const Restrictor &getRestrictor() const
	{
		return *restrictor;
	}
	/**
	 * @brief Set the interpolation operator for interpolating from this level to the finer level.
	 *
	 * @param interpolator the interpolation operator.
	 */
	void setInterpolator(std::shared_ptr<Interpolator> interpolator)
	{
		this->interpolator = interpolator;
	}
	/**
	 * @brief Get the interpolation operator for this level.
	 *
	 * @return Reference to the interpolator.
	 */
	const Interpolator &getInterpolator() const
	{
		return *interpolator;
	}
	/**
	 * @brief Set the operator (matrix) for this level.
	 *
	 * @param op the operator
	 */
	void setOperator(std::shared_ptr<Operator> op)
	{
		this->op = op;
	}
	/**
	 * @brief Get the operator for this level.
	 *
	 * @return Reference to the operator.
	 */
	const Operator &getOperator() const
	{
		return *op;
	}
	/**
	 * @brief Set the smoother for this level.
	 *
	 * @param smoother the smoother
	 */
	void setSmoother(std::shared_ptr<Smoother> smoother)
	{
		this->smoother = smoother;
	}
	/**
	 * @brief Get smoother operator for this level.
	 *
	 * @return Reference to the smoother operator.
	 */
	const Smoother &getSmoother() const
	{
		return *smoother;
	}
	/**
	 * @brief Get DomainCollection for this level.
	 *
	 * @return DomainCollection for this level.
	 */
	const DomainCollection &getDomainCollection() const
	{
		return *dc;
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
