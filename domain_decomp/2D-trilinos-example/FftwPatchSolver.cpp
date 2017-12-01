#include "FftwPatchSolver.h"
using namespace std;

FftwPatchSolver::FftwPatchSolver(DomainSignatureCollection &dsc)
{
	n = dsc.n;
	for (auto &p : dsc.domains) {
		addDomain(p.second);
	}
}
void FftwPatchSolver::addDomain(DomainSignature &d)
{
	if (!initialized) {
		initialized = true;
		f_copy.resize(n * n);
		tmp.resize(n * n);
		sol.resize(n * n);
	}

	if (!plan1.count(d)) {
		fftw_r2r_kind x_transform     = FFTW_RODFT10;
		fftw_r2r_kind x_transform_inv = FFTW_RODFT01;
		fftw_r2r_kind y_transform     = FFTW_RODFT10;
		fftw_r2r_kind y_transform_inv = FFTW_RODFT01;
		if (d.neumann.to_ulong()) {
			if (!d.hasNbr(Side::east) && !d.hasNbr(Side::west)) {
				x_transform     = FFTW_REDFT10;
				x_transform_inv = FFTW_REDFT01;
			} else if (!d.hasNbr(Side::west)) {
				x_transform     = FFTW_REDFT11;
				x_transform_inv = FFTW_REDFT11;
			} else if (!d.hasNbr(Side::east)) {
				x_transform     = FFTW_RODFT11;
				x_transform_inv = FFTW_RODFT11;
			}
			if (!d.hasNbr(Side::north) && !d.hasNbr(Side::south)) {
				y_transform     = FFTW_REDFT10;
				y_transform_inv = FFTW_REDFT01;
			} else if (!d.hasNbr(Side::south)) {
				y_transform     = FFTW_REDFT11;
				y_transform_inv = FFTW_REDFT11;
			} else if (!d.hasNbr(Side::north)) {
				y_transform     = FFTW_RODFT11;
				y_transform_inv = FFTW_RODFT11;
			}
		}

		plan1[d]
		= fftw_plan_r2r_2d(n, n, &f_copy[0], &tmp[0], y_transform, x_transform, FFTW_MEASURE);

		plan2[d]
		= fftw_plan_r2r_2d(n, n, &tmp[0], &sol[0], y_transform_inv, x_transform_inv, FFTW_MEASURE);
	}

	double h_x = d.x_length / n;
	double h_y = d.y_length / n;
	if (!denoms.count(d)) {
		valarray<double> &denom = denoms[d];
		denom.resize(n * n);
		// create denom vector
		if (!d.neumann.to_ulong()) {
			for (int yi = 0; yi < n; yi++) {
				denom[slice(yi * n, n, 1)]
				= -4 / (h_x * h_x) * pow(sin((yi + 1) * M_PI / (2 * n)), 2);
			}
		} else {
			if (!d.hasNbr(Side::north) && !d.hasNbr(Side::south)) {
				for (int yi = 0; yi < n; yi++) {
					denom[slice(yi * n, n, 1)]
					= -4 / (h_x * h_x) * pow(sin(yi * M_PI / (2 * n)), 2);
				}
			} else if (!d.hasNbr(Side::south) || !d.hasNbr(Side::north)) {
				for (int yi = 0; yi < n; yi++) {
					denom[slice(yi * n, n, 1)]
					= -4 / (h_x * h_x) * pow(sin((yi + 0.5) * M_PI / (2 * n)), 2);
				}
			} else {
				for (int yi = 0; yi < n; yi++) {
					denom[slice(yi * n, n, 1)]
					= -4 / (h_x * h_x) * pow(sin((yi + 1) * M_PI / (2 * n)), 2);
				}
			}
		}

		valarray<double> ones(n);
		ones = 1;

		if (d.neumann.none()) {
			for (int xi = 0; xi < n; xi++) {
				denom[slice(xi, n, n)]
				-= 4 / (h_y * h_y) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
			}
		} else {
			if (!d.hasNbr(Side::east) && !d.hasNbr(Side::west)) {
				for (int xi = 0; xi < n; xi++) {
					denom[slice(xi, n, n)]
					-= 4 / (h_y * h_y) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
				}
			} else if (!d.hasNbr(Side::west) || !d.hasNbr(Side::east)) {
				for (int xi = 0; xi < n; xi++) {
					denom[slice(xi, n, n)]
					-= 4 / (h_y * h_y) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
				}
			} else {
				for (int xi = 0; xi < n; xi++) {
					denom[slice(xi, n, n)]
					-= 4 / (h_y * h_y) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
				}
			}
		}
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
void FftwPatchSolver::solve(DomainSignature &d, const vector_type &f, vector_type &u,
                            const vector_type &gamma)
{
	double h_x = d.x_length / n;
	double h_y = d.y_length / n;

	auto f_view     = f.get1dView();
	auto gamma_view = gamma.get1dView();

	int start = d.id * n * n;
	for (int i = 0; i < n * n; i++) {
		f_copy[i] = f_view[start + i];
	}

	if (d.hasNbr(Side::north)) {
		int idx = n * d.index(Side::north);
		for (int i = 0; i < n; i++) {
			f_copy[n * (n - 1) + i] -= 2 / (h_y * h_y) * gamma_view[idx + i];
		}
	}
	if (d.hasNbr(Side::east)) {
		int idx = n * d.index(Side::east);
		for (int i = 0; i < n; i++) {
			f_copy[i * n + (n - 1)] -= 2 / (h_x * h_x) * gamma_view[idx + i];
		}
	}
	if (d.hasNbr(Side::south)) {
		int idx = n * d.index(Side::south);
		for (int i = 0; i < n; i++) {
			f_copy[i] -= 2 / (h_y * h_y) * gamma_view[idx + i];
		}
	}
	if (d.hasNbr(Side::west)) {
		int idx = n * d.index(Side::west);
		for (int i = 0; i < n; i++) {
			f_copy[n * i] -= 2 / (h_x * h_x) * gamma_view[idx + i];
		}
	}

	fftw_execute(plan1[d]);

	tmp /= denoms[d];

	if (d.neumann.all()) {
		tmp[0] = 0;
	}

	fftw_execute(plan2[d]);

	sol /= (4.0 * n * n);

	auto u_view = u.get1dViewNonConst();
	for (int i = 0; i < n * n; i++) {
		u_view[start + i] = sol[i];
	}
}
