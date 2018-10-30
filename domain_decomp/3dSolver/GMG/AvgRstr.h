#ifndef GMGAVGRSTR_H
#define GMGAVGRSTR_H
#include "DomainCollection.h"
#include "InterLevelComm.h"
#include "Restrictor.h"
#include <memory>
namespace GMG
{
/**
 * @brief Restrictor that averages the 8 corresponding fine cells into each coarse cell.
 */
class AvgRstr : public Restrictor
{
	private:
	/**
	 * @brief The coarse DomainCollection that is being restricted to.
	 */
	std::shared_ptr<DomainCollection<3>> coarse_dc;
	/**
	 * @brief The fine DomainCollection that is being restricted from.
	 */
	std::shared_ptr<DomainCollection<3>> fine_dc;
	/**
	 * @brief The communication package for restricting between levels.
	 */
	std::shared_ptr<InterLevelComm> ilc;

	public:
	/**
	 * @brief Create new AvgRstr object.
	 *
	 * @param coarse_dc the DomainColleciton that is being restricted to.
	 * @param fine_dc the DomainCollection that is being restricted from.
	 * @param ilc the communcation package for these two levels.
	 */
	AvgRstr(std::shared_ptr<DomainCollection<3>> coarse_dc,
	        std::shared_ptr<DomainCollection<3>> fine_dc, std::shared_ptr<InterLevelComm> ilc);
	/**
	 * @brief restriction function
	 *
	 * @param coarse the output vector that is restricted to.
	 * @param fine the input vector that is restricted.
	 */
	void restrict(PW<Vec> coarse, PW<Vec> fine) const;
};
} // namespace GMG
#endif
