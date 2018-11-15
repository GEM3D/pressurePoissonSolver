#include "FftwPatchSolver.h"
#include "Utils.h"
using namespace std;
using namespace Utils;
FftwPatchSolver::~FftwPatchSolver()
{
	for (auto p : plan1) {
		fftw_destroy_plan(p.second);
	}
	for (auto p : plan2) {
		fftw_destroy_plan(p.second);
	}
}
void FftwPatchSolver::solve(SchurDomain<3> &d, const Vec f, Vec u, const Vec gamma)
{
	double h_x        = d.domain.lengths[0] / n;
	double h_y        = d.domain.lengths[1] / n;
	double h_z        = d.domain.lengths[2] / n;
	auto   getSpacing = [=](Side<3> s) {
        double retval = 0;
        switch (s.toInt()) {
            case Side<3>::east:
            case Side<3>::west:
                retval = h_x;
                break;
            case Side<3>::south:
            case Side<3>::north:
                retval = h_y;
                break;
            case Side<3>::bottom:
            case Side<3>::top:
                retval = h_z;
        }
        return retval;
	};

	const double *f_view, *gamma_view;
	VecGetArrayRead(f, &f_view);
	VecGetArrayRead(gamma, &gamma_view);

	int start = d.local_index * n * n * n;
	for (int i = 0; i < n * n * n; i++) {
		f_copy[i] = f_view[start + i];
	}

	for (Side<3> s : Side<3>::getValues()) {
		if (d.hasNbr(s)) {
			int    idx = pow(n, 3 - 1) * d.getIfaceLocalIndex(s);
			Slice  sl  = getSlice(&f_copy[0], n, s);
			double h2  = pow(getSpacing(s), 2);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					sl(xi, yi) -= 2.0 / h2 * gamma_view[idx + xi + yi * n];
				}
			}
		}
	}

	fftw_execute(plan1[d]);

	tmp /= denoms[d];

	if (d.neumann.all()) { tmp[0] = 0; }

	fftw_execute(plan2[d]);

	sol /= pow(2.0 * n, 3);

	double *u_view;
	VecGetArray(u, &u_view);
	for (int i = 0; i < n * n * n; i++) {
		u_view[start + i] = sol[i];
	}
	VecRestoreArray(u, &u_view);
	VecRestoreArrayRead(f, &f_view);
	VecRestoreArrayRead(gamma, &gamma_view);
}
