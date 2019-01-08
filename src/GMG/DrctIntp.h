#ifndef GMGDrctIntp_H
#define GMGDrctIntp_H
#include "DomainCollection.h"
#include "InterLevelComm.h"
#include "Interpolator.h"
#include <memory>
namespace GMG
{
/**
 * @brief Simple class that directly places values from coarse vector to fine vector. (This is
 * O(0) accuracy, need to replace with bi-linear interpolation.
 */
template <size_t D> class DrctIntp : public Interpolator
{
	private:
	/**
	 * @brief The coarser set of domains
	 */
	std::shared_ptr<DomainCollection<D>> coarse_dc;
	/**
	 * @brief The finer set of domains
	 */
	std::shared_ptr<DomainCollection<D>> fine_dc;
	/**
	 * @brief The comm package between the levels.
	 */
	std::shared_ptr<InterLevelComm<D>> ilc;

	public:
	/**
	 * @brief Create new DrctIntp object.
	 *
	 * @param coarse_dc the coarser set of domains.
	 * @param fine_dc the finer set of domains.
	 * @param ilc the comm package between the levels.
	 */
	DrctIntp(std::shared_ptr<DomainCollection<D>> coarse_dc,
	         std::shared_ptr<DomainCollection<D>> fine_dc, std::shared_ptr<InterLevelComm<D>> ilc);
	/**
	 * @brief Interpolate from the finer level to the coarser level.
	 *
	 * @param coarse the input vector from the coarser level
	 * @param fine the output vector for the finer level
	 */
	void interpolate(PW<Vec> coarse, PW<Vec> fine) const;
};
template <size_t D>
inline DrctIntp<D>::DrctIntp(std::shared_ptr<DomainCollection<D>> coarse_dc,
                             std::shared_ptr<DomainCollection<D>> fine_dc,
                             std::shared_ptr<InterLevelComm<D>>   ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}
template <size_t D> inline void DrctIntp<D>::interpolate(PW<Vec> coarse, PW<Vec> fine) const
{
	// get vectors
	double *u_fine;
	double *u_coarse;
	PW<Vec> coarse_tmp = ilc->getNewCoarseDistVec();
	// scatter
	PW<VecScatter> scatter = ilc->getScatter();
	VecScatterBegin(scatter, coarse, coarse_tmp, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, coarse, coarse_tmp, INSERT_VALUES, SCATTER_FORWARD);

	VecGetArray(fine, &u_fine);
	VecGetArray(coarse_tmp, &u_coarse);
	for (auto p : ilc->getFineDomains()) {
		Domain<D> &d          = *p.d;
		int        n          = d.n;
		int        coarse_idx = p.local_index * pow(n, D);
		int        fine_idx   = d.id_local * pow(n, D);

		Orthant<D> orth = d.oct_on_parent;
		if (d.id != d.parent_id) {
			std::array<int, D> strides;
			std::array<int, D> starts;
			for (size_t i = 0; i < D; i++) {
				strides[i] = pow(n, i);
				starts[i]  = orth.isOnSide(2 * i) ? 0 : n;
			}
			for (int i = 0; i < (int) pow(n, D); i++) {
				std::array<int, D> coord;
				int                idx = 0;
				for (size_t x = 0; x < D; x++) {
					coord[x] = (i / strides[x]) % n;
					idx += (coord[x] + starts[x]) / 2 * strides[x];
				}
				u_fine[fine_idx + i] += u_coarse[coarse_idx + idx];
			}
		} else {
			for (int i = 0; i < pow(n, D); i++) {
				u_fine[fine_idx + i] += u_coarse[coarse_idx + i];
			}
		}
	}
	VecRestoreArray(fine, &u_fine);
	VecRestoreArray(coarse_tmp, &u_coarse);
}
} // namespace GMG
#endif
