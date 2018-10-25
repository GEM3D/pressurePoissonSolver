#include "FftwPatchSolver.h"
#include "Utils.h"
using namespace std;
using namespace Utils;
inline int index(const int &n, const int &xi, const int &yi, const int &zi)
{
	return xi + yi * n + zi * n * n;
}
FftwPatchSolver::FftwPatchSolver(DomainCollection &dc, double lambda)
{
	n            = dc.getN();
	this->lambda = lambda;
}
void FftwPatchSolver::addDomain(SchurDomain<3> &d)
{
	if (!initialized) {
		initialized = true;
		f_copy.resize(n * n * n);
		tmp.resize(n * n * n);
		sol.resize(n * n * n);
	}

	if (!plan1.count(d)) {
		fftw_r2r_kind x_transform     = FFTW_RODFT10;
		fftw_r2r_kind x_transform_inv = FFTW_RODFT01;
		fftw_r2r_kind y_transform     = FFTW_RODFT10;
		fftw_r2r_kind y_transform_inv = FFTW_RODFT01;
		fftw_r2r_kind z_transform     = FFTW_RODFT10;
		fftw_r2r_kind z_transform_inv = FFTW_RODFT01;
		if (d.neumann.to_ulong()) {
			// x direction
			if (d.isNeumann(Side::east) && d.isNeumann(Side::west)) {
				x_transform     = FFTW_REDFT10;
				x_transform_inv = FFTW_REDFT01;
			} else if (d.isNeumann(Side::west)) {
				x_transform     = FFTW_REDFT11;
				x_transform_inv = FFTW_REDFT11;
			} else if (d.isNeumann(Side::east)) {
				x_transform     = FFTW_RODFT11;
				x_transform_inv = FFTW_RODFT11;
			}
			// y direction
			if (d.isNeumann(Side::north) && d.isNeumann(Side::south)) {
				y_transform     = FFTW_REDFT10;
				y_transform_inv = FFTW_REDFT01;
			} else if (d.isNeumann(Side::south)) {
				y_transform     = FFTW_REDFT11;
				y_transform_inv = FFTW_REDFT11;
			} else if (d.isNeumann(Side::north)) {
				y_transform     = FFTW_RODFT11;
				y_transform_inv = FFTW_RODFT11;
			}
			// z direction
			if (d.isNeumann(Side::bottom) && d.isNeumann(Side::top)) {
				z_transform     = FFTW_REDFT10;
				z_transform_inv = FFTW_REDFT01;
			} else if (d.isNeumann(Side::bottom)) {
				z_transform     = FFTW_REDFT11;
				z_transform_inv = FFTW_REDFT11;
			} else if (d.isNeumann(Side::top)) {
				z_transform     = FFTW_RODFT11;
				z_transform_inv = FFTW_RODFT11;
			}
		}

		plan1[d] = fftw_plan_r2r_3d(n, n, n, &f_copy[0], &tmp[0], z_transform, y_transform,
		                            x_transform, FFTW_MEASURE | FFTW_DESTROY_INPUT);

		plan2[d] = fftw_plan_r2r_3d(n, n, n, &tmp[0], &sol[0], z_transform_inv, y_transform_inv,
		                            x_transform_inv, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	}

	double h_x = d.domain.lengths[0] / n;
	double h_y = d.domain.lengths[1] / n;
	double h_z = d.domain.lengths[2] / n;
	if (!denoms.count(d)) {
		valarray<double> &denom = denoms[d];
		denom.resize(n * n * n);
		// create denom vector
		// z direction
		if (d.isNeumann(Side::bottom) && d.isNeumann(Side::top)) {
			for (int zi = 0; zi < n; zi++) {
				denom[slice(zi * n * n, n * n, 1)]
				= -4 / (h_z * h_z) * pow(sin(zi * M_PI / (2 * n)), 2);
			}
		} else if (d.isNeumann(Side::bottom) || d.isNeumann(Side::top)) {
			for (int zi = 0; zi < n; zi++) {
				denom[slice(zi * n * n, n * n, 1)]
				= -4 / (h_z * h_z) * pow(sin((zi + 0.5) * M_PI / (2 * n)), 2);
			}
		} else {
			for (int zi = 0; zi < n; zi++) {
				denom[slice(zi * n * n, n * n, 1)]
				= -4 / (h_z * h_z) * pow(sin((zi + 1) * M_PI / (2 * n)), 2);
			}
		}

		valarray<double> ones(n * n);
		ones = 1;

		// y direction
		valarray<size_t> sizes = {(size_t) n, (size_t) n};
		valarray<size_t> strides(2);
		strides[0] = n * n;
		strides[1] = 1;
		if (d.isNeumann(Side::south) && d.isNeumann(Side::north)) {
			for (int yi = 0; yi < n; yi++) {
				denom[gslice(yi * n, sizes, strides)]
				-= 4 / (h_y * h_y) * pow(sin(yi * M_PI / (2 * n)), 2) * ones;
			}
		} else if (d.isNeumann(Side::south) || d.isNeumann(Side::north)) {
			for (int yi = 0; yi < n; yi++) {
				denom[gslice(yi * n, sizes, strides)]
				-= 4 / (h_y * h_y) * pow(sin((yi + 0.5) * M_PI / (2 * n)), 2) * ones;
			}
		} else {
			for (int yi = 0; yi < n; yi++) {
				denom[gslice(yi * n, sizes, strides)]
				-= 4 / (h_y * h_y) * pow(sin((yi + 1) * M_PI / (2 * n)), 2) * ones;
			}
		}

		// x direction
		strides[0] = n * n;
		strides[1] = n;
		if (d.isNeumann(Side::west) && d.isNeumann(Side::east)) {
			for (int xi = 0; xi < n; xi++) {
				denom[gslice(xi, sizes, strides)]
				-= 4 / (h_x * h_x) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
			}
		} else if (d.isNeumann(Side::west) || d.isNeumann(Side::east)) {
			for (int xi = 0; xi < n; xi++) {
				denom[gslice(xi, sizes, strides)]
				-= 4 / (h_x * h_x) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
			}
		} else {
			for (int xi = 0; xi < n; xi++) {
				denom[gslice(xi, sizes, strides)]
				-= 4 / (h_x * h_x) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
			}
		}
		denom += lambda;
	}
}
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
	auto   getSpacing = [=](Side s) {
        double retval = 0;
        switch (s.toInt()) {
            case Side::east:
            case Side::west:
                retval = h_x;
                break;
            case Side::south:
            case Side::north:
                retval = h_y;
                break;
            case Side::bottom:
            case Side::top:
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

	for (Side s : Side::getValues()) {
		if (d.hasNbr(s)) {
			int    idx = n * n * d.getIfaceLocalIndex(s);
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

	sol /= (8.0 * n * n * n);

	double *u_view;
	VecGetArray(u, &u_view);
	for (int i = 0; i < n * n * n; i++) {
		u_view[start + i] = sol[i];
	}
	VecRestoreArray(u, &u_view);
	VecRestoreArrayRead(f, &f_view);
	VecRestoreArrayRead(gamma, &gamma_view);
}
