#ifndef InterLevelComm_H
#define InterLevelComm_H
#include "DomainCollection.h"
namespace GMG
{
/**
 * @brief Structure that wraps some extra meta-data around a Domain object.
 */
struct ILCFineToCoarseMetadata {
	/**
	 * @brief The domain that this meta-data corresponds to.
	 */
	std::shared_ptr<Domain<3>> d;
	/**
	 * @brief the local block-index of the parent domain in the scattered coarse vector.
	 */
	int local_index;
	/**
	 * @brief the global block-index of the parent domain in the scattered coarse vector.
	 */
	int global_index;
	/**
	 * @brief less than operator so that this struct can be placed in set container.
	 */
	bool operator<(const ILCFineToCoarseMetadata &other) const
	{
		return *d < *other.d;
	}
};
/**
 * @brief Creates a mapping from fine to coarse levels.
 */
class InterLevelComm
{
	private:
	/**
	 * @brief The number of cells in each direction on a patch
	 */
	int n;
	/**
	 * @brief the number of elements in the local distributed vector.
	 */
	int local_vec_size;
	/**
	 * @brief a set of domains on the finer level that contains extra meta-data on the mapping from
	 * finer level vectors to coarser level vectors.
	 */
	std::set<ILCFineToCoarseMetadata> coarse_domains;
	/**
	 * @brief The PETSc VecScatter object that scatters the values of the coarse vector to each
	 * process that needs them for interpolation / restriction.
	 */
	PW<VecScatter> scatter;

	public:
	/**
	 * @brief Create a new InterLevelComm object.
	 *
	 * @param coarse_dc the coarser DomainCollection.
	 * @param fine_dc the finer DomainCollection.
	 */
	InterLevelComm(std::shared_ptr<DomainCollection> coarse_dc,
	               std::shared_ptr<DomainCollection> fine_dc);

	/**
	 * @brief Allocate a new vector for which coarse values values will be scattered into.
	 *
	 * @return the newly allocated vector.
	 */
	PW_explicit<Vec> getNewCoarseDistVec();

	/**
	 * @brief Get a PETSc VecScatter object that will scatter values from the coarse vector to the
	 * processes that need them for interpolation / restriction.
	 *
	 * @return the VecScatter object
	 */
	PW_explicit<VecScatter> getScatter();

	/**
	 * @brief get a set of domains on the finer level that contains meta-data for how the values in
	 * the fine vector map to the values in the coarse vector.
	 *
	 * @return set of fine domains
	 */
	std::set<ILCFineToCoarseMetadata> getFineDomains()
	{
		return coarse_domains;
	}
};
} // namespace GMG
#endif
