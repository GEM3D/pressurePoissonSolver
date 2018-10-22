#include "AvgRstr.h"
#include <functional>
using namespace std;
using namespace GMG;
AvgRstr::AvgRstr(shared_ptr<DomainCollection> coarse_dc, shared_ptr<DomainCollection> fine_dc,
                 shared_ptr<InterLevelComm> ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}
void AvgRstr::restrict(PW<Vec> coarse, PW<Vec> fine) const
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
	for (ILCFineToCoarseMetadata data : ilc->getFineDomains()) {
		Domain<3> &                             d          = *data.d;
		int                                     n          = d.n;
		int                                     coarse_idx = data.local_index * n * n * n;
		const function<double &(int, int, int)> f_to_c[8]  = {
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi) / 2) + ((yi) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi + n) / 2) + ((yi) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi) / 2) + ((yi + n) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi + n) / 2) + ((yi + n) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi) / 2) + ((yi) / 2) * n + ((zi + n) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi + n) / 2) + ((yi) / 2) * n + ((zi + n) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi) / 2) + ((yi + n) / 2) * n + ((zi + n) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return f_coarse[coarse_idx + ((xi + n) / 2) + ((yi + n) / 2) * n
                            + ((zi + n) / 2) * n * n];
        }};

		int fine_idx = d.id_local * n * n * n;
		int oct      = d.oct_on_parent;
		if (d.id != d.parent_id) {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						f_to_c[oct](xi, yi, zi) += r_fine[fine_idx + xi + yi * n + zi * n * n] / 8;
					}
				}
			}
		} else {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						f_coarse[coarse_idx + xi + yi * n + zi * n * n]
						+= r_fine[fine_idx + xi + yi * n + zi * n * n];
					}
				}
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
