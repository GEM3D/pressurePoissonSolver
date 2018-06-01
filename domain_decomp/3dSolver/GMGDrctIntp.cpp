#include "GMGDrctIntp.h"
#include <functional>
using namespace std;
GMGDrctIntp::GMGDrctIntp(DomainCollection &coarse_dc, DomainCollection &fine_dc)
{
	this->coarse_dc = coarse_dc;
	this->fine_dc   = fine_dc;
}
void GMGDrctIntp::interpolate(PW<Vec> coarse, PW<Vec> fine)
{
	// get vectors
	double *u_fine;
	double *u_coarse;
	VecGetArray(fine, &u_fine);
	VecGetArray(coarse, &u_coarse);
	for (auto p : coarse_dc.domains) {
		Domain &                                d          = p.second;
		int                                     n          = d.n;
		int                                     coarse_idx = d.id_local * n * n * n;
		const function<double &(int, int, int)> f_to_c[8]  = {
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi) / 2) + ((yi) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi + n) / 2) + ((yi) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi) / 2) + ((yi + n) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi + n) / 2) + ((yi + n) / 2) * n + ((zi) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi) / 2) + ((yi) / 2) * n + ((zi + n) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi + n) / 2) + ((yi) / 2) * n + ((zi + n) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi) / 2) + ((yi + n) / 2) * n + ((zi + n) / 2) * n * n];
        },
        [&](int xi, int yi, int zi) -> double & {
            return u_coarse[coarse_idx + ((xi + n) / 2) + ((yi + n) / 2) * n
                            + ((zi + n) / 2) * n * n];
        }};
		for (int oct = 0; oct < 8; oct++) {
			int fine_idx = fine_dc.domains.at(d.child_id[oct]).id_local * n * n * n;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n] += f_to_c[oct](xi, yi, zi);
					}
				}
			}
		}
	}
	VecRestoreArray(fine, &u_fine);
	VecRestoreArray(coarse, &u_coarse);
}
