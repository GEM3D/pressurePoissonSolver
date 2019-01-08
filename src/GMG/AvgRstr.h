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
template <size_t D> class AvgRstr : public Restrictor
{
	private:
	/**
	 * @brief The coarse DomainCollection that is being restricted to.
	 */
	std::shared_ptr<DomainCollection<D>> coarse_dc;
	/**
	 * @brief The fine DomainCollection that is being restricted from.
	 */
	std::shared_ptr<DomainCollection<D>> fine_dc;
	/**
	 * @brief The communication package for restricting between levels.
	 */
	std::shared_ptr<InterLevelComm<D>> ilc;

	public:
	/**
	 * @brief Create new AvgRstr object.
	 *
	 * @param coarse_dc the DomainColleciton that is being restricted to.
	 * @param fine_dc the DomainCollection that is being restricted from.
	 * @param ilc the communcation package for these two levels.
	 */
	AvgRstr(std::shared_ptr<DomainCollection<D>> coarse_dc,
	        std::shared_ptr<DomainCollection<D>> fine_dc, std::shared_ptr<InterLevelComm<D>> ilc);
	/**
	 * @brief restriction function
	 *
	 * @param coarse the output vector that is restricted to.
	 * @param fine the input vector that is restricted.
	 */
	void restrict(PW<Vec> coarse, PW<Vec> fine) const;
};

template <size_t D>
inline AvgRstr<D>::AvgRstr(std::shared_ptr<DomainCollection<D>> coarse_dc,
                           std::shared_ptr<DomainCollection<D>> fine_dc,
                           std::shared_ptr<InterLevelComm<D>>   ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}
template <size_t D> inline void AvgRstr<D>::restrict(PW<Vec> coarse, PW<Vec> fine) const
{
	// get vectors
	VecSet(coarse, 0);
	double *r_fine;
	double *f_coarse;
	// store in tmp vector for fine level
	PW<Vec> coarse_tmp = ilc->getNewCoarseDistVec();
	VecSet(coarse_tmp, 0);
	VecGetArray(fine, &r_fine);
	VecGetArray(coarse_tmp, &f_coarse);
	for (ILCFineToCoarseMetadata<D> data : ilc->getFineDomains()) {
		Domain<D> &d          = *data.d;
		int        n          = d.n;
		int        coarse_idx = data.local_index * pow(n, D);
		int        fine_idx   = d.id_local * pow(n, D);
		Orthant<D> orth       = d.oct_on_parent;
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
				f_coarse[coarse_idx + idx] += r_fine[fine_idx + i] / (1 << D);
			}
		} else {
			for (int i = 0; i < pow(n, D); i++) {
				f_coarse[coarse_idx + i] += r_fine[fine_idx + i];
			}
		}
	}
	VecRestoreArray(fine, &r_fine);
	VecRestoreArray(coarse_tmp, &f_coarse);
	// scatter
	PW<VecScatter> scatter = ilc->getScatter();
	VecScatterBegin(scatter, coarse_tmp, coarse, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, coarse_tmp, coarse, ADD_VALUES, SCATTER_REVERSE);
}
} // namespace GMG
#endif
