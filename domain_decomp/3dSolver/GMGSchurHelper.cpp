#include "GMGSchurHelper.h"
using namespace std;
GMGSchurHelper::GMGSchurHelper(int n, OctTree t, DomainCollection &dc, SchurHelper &sh)
{
	num_levels = t.num_levels;
	top_level  = t.num_levels - 1;
	levels.resize(num_levels);
	shs.resize(num_levels);
	u_vectors.resize(num_levels);
	f_vectors.resize(num_levels);
	r_vectors.resize(num_levels);
	b_vectors.resize(num_levels);
	g_vectors.resize(num_levels);
	jacobis.resize(num_levels);
	levels[top_level]    = dc;
	shs[top_level]       = sh;
	u_vectors[top_level] = dc.getNewDomainVec();
	f_vectors[top_level] = dc.getNewDomainVec();
	r_vectors[top_level] = dc.getNewDomainVec();
	jacobis[top_level]   = shs[top_level].formPBMatrix()->getBlockJacobiSmoother();
	for (int i = 0; i < top_level; i++) {
		levels[i]   = DomainCollection(t, i + 1,n);
		shs[i]       = SchurHelper(levels[i], sh.getSolver(), sh.getOp(), sh.getInterpolator());
		f_vectors[i] = levels[i].getNewDomainVec();
		u_vectors[i] = levels[i].getNewDomainVec();
		r_vectors[i] = levels[i].getNewDomainVec();
		b_vectors[i] = shs[i].getNewSchurVec();
		g_vectors[i] = shs[i].getNewSchurVec();
		if (i != 0) { jacobis[i] = shs[i].formPBMatrix()->getBlockJacobiSmoother(); }
	}
}
void GMGSchurHelper::restrictForLevel(int level)
{
	const DomainCollection &coarse_dc = levels[level];
	const DomainCollection &fine_dc   = levels[level + 1];
	// get vectors
	VecScale(f_vectors[level], 0);
	double *r_fine;
	double *f_coarse;
	VecGetArray(r_vectors[level + 1], &r_fine);
	VecGetArray(f_vectors[level], &f_coarse);
	for (auto p : coarse_dc.domains) {
		Domain &                                d          = *p.second;
		int                                     n          = d.n;
		int                                     coarse_idx = d.id_local * n * n * n;
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

		for (int oct = 0; oct < 8; oct++) {
			int fine_idx = fine_dc.domains.at(d.child_id[oct])->id_local * n * n * n;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						f_to_c[oct](xi, yi, zi) += r_fine[fine_idx + xi + yi * n + zi * n * n] / 8;
					}
				}
			}
		}
	}
	VecRestoreArray(r_vectors[level + 1], &r_fine);
	VecRestoreArray(f_vectors[level], &f_coarse);
}
void GMGSchurHelper::prolongateFromLevel(int level)
{
	const DomainCollection &coarse_dc = levels[level];
	const DomainCollection &fine_dc   = levels[level + 1];
	// get vectors
	double *u_fine;
	double *u_coarse;
	VecGetArray(u_vectors[level + 1], &u_fine);
	VecGetArray(u_vectors[level], &u_coarse);
	for (auto p : coarse_dc.domains) {
		Domain &                                d          = *p.second;
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
			int fine_idx = fine_dc.domains.at(d.child_id[oct])->id_local * n * n * n;
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					for (int xi = 0; xi < n; xi++) {
						u_fine[fine_idx + xi + yi * n + zi * n * n] += f_to_c[oct](xi, yi, zi);
					}
				}
			}
		}
	}
	VecRestoreArray(u_vectors[level + 1], &u_fine);
	VecRestoreArray(u_vectors[level], &u_coarse);
}
void GMGSchurHelper::apply(Vec b, Vec gamma)
{
	// finest level
	jacobis[top_level].apply(gamma, b);
	jacobis[top_level].apply(gamma, b);
	shs[top_level].apply(u_vectors[top_level], r_vectors[top_level]);
	VecAYPX(r_vectors[top_level], -1, f_vectors[top_level]);
	// down-cycle
	for (int i = top_level - 1; i >= 1; i--) {
		restrictForLevel(i);
		VecScale(u_vectors[i], 0);
		VecScale(b_vectors[i], 0);
		VecScale(g_vectors[i], 0);
		shs[i].solveWithInterface(f_vectors[i], u_vectors[i], g_vectors[i], b_vectors[i]);
		VecScale(b_vectors[i], -1.0);
		jacobis[i].apply(g_vectors[i], b_vectors[i]);
		shs[i].apply(u_vectors[i], r_vectors[i]);
		VecAYPX(r_vectors[i], -1, f_vectors[i]);
	}
	// coarse level
	VecScale(u_vectors[0], 0);
	VecScale(u_vectors[0], 0);
	VecScale(b_vectors[0], 0);
	VecScale(g_vectors[0], 0);
	restrictForLevel(0);
	shs[0].solveWithInterface(f_vectors[0], u_vectors[0], g_vectors[0], b_vectors[0]);
	prolongateFromLevel(0);
	// up cycle
	for (int i = 1; i <= num_levels - 2; i++) {
		//	shs[i].interpolateToInterface(f_vectors[i], u_vectors[i], g_vectors[i]);
		// shs[i].solveWithSolution(f_vectors[i], u_vectors[i]);
		// shs[i].solveWithInterface(f_vectors[i], u_vectors[i], g_vectors[i], b_vectors[i]);
		prolongateFromLevel(i);
	}
	// finest level
	// shs[top_level].interpolateToInterface(f_vectors[top_level], u_vectors[top_level],
	//                                    gamma);
}
