#include "GMGHelper.h"
using namespace std;
GMGHelper::GMGHelper(int n, OctTree t, std::shared_ptr<PatchSolver> solver,
                     std::shared_ptr<PatchOperator> op, std::shared_ptr<Interpolator> interpolator)
{
	num_levels = t.num_levels;
	levels.resize(num_levels);
	smoothers.resize(num_levels);
	shs.resize(num_levels);
	u_vectors.resize(num_levels);
	f_vectors.resize(num_levels);
	r_vectors.resize(num_levels);
	for (int i = 0; i < num_levels; i++) {
		levels[i]   = DomainCollection(t, i + 1);
		levels[i].n = n;
		for (auto &p : levels[i].domains) {
			p.second.n = n;
			solver->addDomain(p.second);
		}
		shs[i]       = SchurHelper(levels[i], solver, op, interpolator);
		smoothers[i] = SchwarzPrec(&shs[i], &levels[i]);
		f_vectors[i] = levels[i].getNewDomainVec();
		u_vectors[i] = levels[i].getNewDomainVec();
		r_vectors[i] = levels[i].getNewDomainVec();
	}
}
void GMGHelper::restrictForLevel(int level)
{
	DomainCollection &coarse_dc = levels[level];
	DomainCollection &fine_dc   = levels[level + 1];
	// get vectors
	VecScale(f_vectors[level], 0);
	double *r_fine;
	double *f_coarse;
	VecGetArray(r_vectors[level + 1], &r_fine);
	VecGetArray(f_vectors[level], &f_coarse);
	for (auto &p : coarse_dc.domains) {
		Domain &d          = p.second;
		int     n          = d.n;
		int     coarse_idx = d.id_local * n * n * n;
		const function<double &(int, int, int)> f_to_c[8]
		= {[&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi) / 2 + (yi) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi + n) / 2 + (yi) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi) / 2 + (yi + n) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi + n) / 2 + (yi + n) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi) / 2 + (yi) / 2 * n + (zi + n) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi + n) / 2 + (yi) / 2 * n + (zi + n) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi) / 2 + (yi + n) / 2 * n + (zi + n) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return f_coarse[coarse_idx + (xi + n) / 2 + (yi + n) / 2 * n + (zi + n) / 2 * n * n];
		   }};
		for (int oct = 0; oct < 8; oct++) {
			int fine_idx = fine_dc.domains[d.child_id[oct]].id_local * n * n * n;
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
void GMGHelper::prolongateFromLevel(int level)
{
	DomainCollection &coarse_dc = levels[level];
	DomainCollection &fine_dc   = levels[level + 1];
	// get vectors
	double *u_fine;
	double *u_coarse;
	VecGetArray(u_vectors[level + 1], &u_fine);
	VecGetArray(u_vectors[level], &u_coarse);
	for (auto &p : coarse_dc.domains) {
		Domain &d          = p.second;
		int     n          = d.n;
		int     coarse_idx = d.id_local * n * n * n;
		const function<double &(int, int, int)> f_to_c[8]
		= {[&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi) / 2 + (yi) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi + n) / 2 + (yi) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi) / 2 + (yi + n) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi + n) / 2 + (yi + n) / 2 * n + (zi) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi) / 2 + (yi) / 2 * n + (zi + n) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi + n) / 2 + (yi) / 2 * n + (zi + n) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi) / 2 + (yi + n) / 2 * n + (zi + n) / 2 * n * n];
		   },
		   [&](int xi, int yi, int zi) -> double & {
			   return u_coarse[coarse_idx + (xi + n) / 2 + (yi + n) / 2 * n + (zi + n) / 2 * n * n];
		   }};
		for (int oct = 0; oct < 8; oct++) {
			int fine_idx = fine_dc.domains[d.child_id[oct]].id_local * n * n * n;
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
void GMGHelper::apply(Vec f, Vec u)
{
	f_vectors[num_levels - 1] = f;
	u_vectors[num_levels - 1] = u;

	// finest level
	smoothers[num_levels - 1].apply(f, u);
	shs[num_levels - 1].apply(u_vectors[num_levels - 1], r_vectors[num_levels - 1]);
	VecAYPX(r_vectors[num_levels - 1], -1, f_vectors[num_levels - 1]);
	// down-cycle
	for (int i = num_levels - 2; i > 0; i--) {
		restrictForLevel(i);
		VecScale(u_vectors[i], 0);
		smoothers[i].apply(f_vectors[i], u_vectors[i]);
		shs[i].apply(u_vectors[i], r_vectors[i]);
		VecAYPX(r_vectors[i], -1, f_vectors[i]);
	}
	// coarse level
	restrictForLevel(0);
	VecScale(u_vectors[0], 0);
	smoothers[0].apply(f_vectors[0], u_vectors[0]);
	prolongateFromLevel(0);
	// up cycle
	for (int i = 1; i < num_levels - 1; i++) {
		smoothers[i].apply(f_vectors[i], u_vectors[i]);
		prolongateFromLevel(i);
	}
	// finest level
	smoothers[num_levels - 1].apply(f_vectors[num_levels - 1], u_vectors[num_levels - 1]);
}
