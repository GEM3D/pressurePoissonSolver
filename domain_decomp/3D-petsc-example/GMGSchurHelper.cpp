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
	levels[top_level]    = dc;
	shs[top_level]       = sh;
	r_vectors[top_level] = dc.getNewDomainVec();
	for (int i = 0; i < top_level; i++) {
		levels[i]   = DomainCollection(t, i + 1);
		levels[i].n = n;
		for (auto &p : levels[i].domains) {
			p.second.n = n;
		}
		shs[i]       = SchurHelper(levels[i], sh.getSolver(), sh.getOp(), sh.getInterpolator());
		f_vectors[i] = levels[i].getNewDomainVec();
		u_vectors[i] = levels[i].getNewDomainVec();
		r_vectors[i] = levels[i].getNewDomainVec();
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
		Domain &d          = p.second;
		int     n          = d.n;
		int     coarse_idx = d.id_local * n * n * n;
		const function<double &(int, int, int)> f_to_c[8] = {
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
			int fine_idx = fine_dc.domains.at(d.child_id[oct]).id_local * n * n * n;
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
		Domain &d          = p.second;
		int     n          = d.n;
		int     coarse_idx = d.id_local * n * n * n;
		const function<double &(int, int, int)> f_to_c[8] = {
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
	VecRestoreArray(u_vectors[level + 1], &u_fine);
	VecRestoreArray(u_vectors[level], &u_coarse);
}
void GMGSchurHelper::apply(Vec b, Vec gamma)
{
	// finest level
	//shs[top_level].apply(u, r_vectors[top_level]);
	//VecAYPX(r_vectors[top_level], -1, f);
	// down-cycle
	for (int i = top_level - 1; i >= 1; i--) {
		restrictForLevel(i);
		VecScale(u_vectors[i], 0);
		shs[i].solveWithSolution(f_vectors[i], u_vectors[i]);
		shs[i].apply(u_vectors[i], r_vectors[i]);
		VecAYPX(r_vectors[i], -1, f_vectors[i]);
	}
	// coarse level
	restrictForLevel(0);
	VecScale(u_vectors[0], 0);
	shs[0].solveWithSolution(f_vectors[0], u_vectors[0]);
	prolongateFromLevel(0);
	// up cycle
	for (int i = 1; i <= num_levels - 2; i++) {
		shs[i].solveWithSolution(f_vectors[i], u_vectors[i]);
		prolongateFromLevel(i);
	}
	// finest level
}
