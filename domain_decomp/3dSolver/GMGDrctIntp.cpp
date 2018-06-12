#include "GMGDrctIntp.h"
#include <functional>
using namespace std;
GMGDrctIntp::GMGDrctIntp(shared_ptr<DomainCollection> coarse_dc,
                         shared_ptr<DomainCollection> fine_dc, shared_ptr<InterLevelComm> ilc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
	this->ilc       = ilc;
}
void GMGDrctIntp::interpolate(PW<Vec> coarse, PW<Vec> fine)
{
	// get vectors
	double *u_fine;
	double *u_coarse;
	PW<Vec> coarse_tmp = ilc->getNewCoarseDistVec();
	// scatter
	PW<VecScatter> scatter = ilc->getScatter();
	VecScatterBegin(scatter, coarse, coarse_tmp, ADD_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, coarse, coarse_tmp, ADD_VALUES, SCATTER_FORWARD);

	VecGetArray(fine, &u_fine);
	VecGetArray(coarse_tmp, &u_coarse);
	for (auto p : ilc->getFineDomains()) {
		Domain &d          = *p.d;
		int     n          = d.n;
		int     coarse_idx = p.local_index * n * n * n;

		const function<void(double *, int, double *, int)> f_to_c[8] = {
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi) / 2) + ((yi) / 2) * n + ((zi) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi + n) / 2) + ((yi) / 2) * n
						            + ((zi) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi) / 2) + ((yi + n) / 2) * n
						            + ((zi) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi + n) / 2) + ((yi + n) / 2) * n
						            + ((zi) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi) / 2) + ((yi) / 2) * n
						            + ((zi + n) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi + n) / 2) + ((yi) / 2) * n
						            + ((zi + n) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi) / 2) + ((yi + n) / 2) * n
						            + ((zi + n) / 2) * n * n];
					}
				}
			}
		},
		[&](double *u_fine, int fine_idx, double *u_coarse, int coarse_idx) -> void {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n]
						+= u_coarse[coarse_idx + ((xi + n) / 2) + ((yi + n) / 2) * n
						            + ((zi + n) / 2) * n * n];
					}
				}
			}
		}};
		int oct      = d.oct_on_parent;
		int fine_idx = d.id_local * n * n * n;
		f_to_c[oct](u_fine, fine_idx, u_coarse, coarse_idx);
	}
	VecRestoreArray(fine, &u_fine);
	VecRestoreArray(coarse_tmp, &u_coarse);
}
