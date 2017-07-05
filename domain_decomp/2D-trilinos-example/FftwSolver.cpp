#include "FftwSolver.h"
using namespace std;
FftwSolver::FftwSolver(Domain *dom)
{
	this->d = dom;
	f_copy  = valarray<double>(d->n * d->n);
	tmp     = valarray<double>(d->n * d->n);
	u       = valarray<double>(d->n * d->n);
	denom   = valarray<double>(d->n * d->n);

	fftw_r2r_kind x_transform     = FFTW_RODFT10;
	fftw_r2r_kind x_transform_inv = FFTW_RODFT01;
	fftw_r2r_kind y_transform     = FFTW_RODFT10;
	fftw_r2r_kind y_transform_inv = FFTW_RODFT01;
	if (d->neumann) {
		if (!d->hasNbr(Side::east) && !d->hasNbr(Side::west)) {
			x_transform     = FFTW_REDFT10;
			x_transform_inv = FFTW_REDFT01;
		} else if (!d->hasNbr(Side::west)) {
			x_transform     = FFTW_REDFT11;
			x_transform_inv = FFTW_REDFT11;
		} else if (!d->hasNbr(Side::east)) {
			x_transform     = FFTW_RODFT11;
			x_transform_inv = FFTW_RODFT11;
		}
		if (!d->hasNbr(Side::north) && !d->hasNbr(Side::south)) {
			y_transform     = FFTW_REDFT10;
			y_transform_inv = FFTW_REDFT01;
		} else if (!d->hasNbr(Side::south)) {
			y_transform     = FFTW_REDFT11;
			y_transform_inv = FFTW_REDFT11;
		} else if (!d->hasNbr(Side::north)) {
			y_transform     = FFTW_RODFT11;
			y_transform_inv = FFTW_RODFT11;
		}
	}

	plan1
	= fftw_plan_r2r_2d(d->n, d->n, &f_copy[0], &tmp[0], y_transform, x_transform, FFTW_MEASURE);

	plan2
	= fftw_plan_r2r_2d(d->n, d->n, &tmp[0], &u[0], y_transform_inv, x_transform_inv, FFTW_MEASURE);

	// create denom vector
	if (!d->neumann) {
		for (int yi = 0; yi < d->n; yi++) {
			denom[slice(yi * d->n, d->n, 1)]
			= -4 / (d->h_x * d->h_x) * pow(sin((yi + 1) * M_PI / (2 * d->n)), 2);
		}
	} else {
		if (!d->hasNbr(Side::north) && !d->hasNbr(Side::south)) {
			for (int yi = 0; yi < d->n; yi++) {
				denom[slice(yi * d->n, d->n, 1)]
				= -4 / (d->h_x * d->h_x) * pow(sin(yi * M_PI / (2 * d->n)), 2);
			}
		} else if (!d->hasNbr(Side::south) || !d->hasNbr(Side::north)) {
			for (int yi = 0; yi < d->n; yi++) {
				denom[slice(yi * d->n, d->n, 1)]
				= -4 / (d->h_x * d->h_x) * pow(sin((yi + 0.5) * M_PI / (2 * d->n)), 2);
			}
		} else {
			for (int yi = 0; yi < d->n; yi++) {
				denom[slice(yi * d->n, d->n, 1)]
				= -4 / (d->h_x * d->h_x) * pow(sin((yi + 1) * M_PI / (2 * d->n)), 2);
			}
		}
	}

	valarray<double> ones(d->n);
	ones = 1;

	if (!d->neumann) {
		for (int xi = 0; xi < d->n; xi++) {
			denom[slice(xi, d->n, d->n)]
			-= 4 / (d->h_y * d->h_y) * pow(sin((xi + 1) * M_PI / (2 * d->n)), 2) * ones;
		}
	} else {
		if (!d->hasNbr(Side::east) && !d->hasNbr(Side::west)) {
			for (int xi = 0; xi < d->n; xi++) {
				denom[slice(xi, d->n, d->n)]
				-= 4 / (d->h_y * d->h_y) * pow(sin(xi * M_PI / (2 * d->n)), 2) * ones;
			}
		} else if (!d->hasNbr(Side::west) || !d->hasNbr(Side::east)) {
			for (int xi = 0; xi < d->n; xi++) {
				denom[slice(xi, d->n, d->n)]
				-= 4 / (d->h_y * d->h_y) * pow(sin((xi + 0.5) * M_PI / (2 * d->n)), 2) * ones;
			}
		} else {
			for (int xi = 0; xi < d->n; xi++) {
				denom[slice(xi, d->n, d->n)]
				-= 4 / (d->h_y * d->h_y) * pow(sin((xi + 1) * M_PI / (2 * d->n)), 2) * ones;
			}
		}
	}
}
FftwSolver::~FftwSolver()
{
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
}
void FftwSolver::solve()
{
	f_copy = d->f;
	if (!d->hasNbr(Side::north) && d->neumann) {
		f_copy[slice(d->n * (d->n - 1), d->n, 1)] -= 1 / d->h_y * d->boundary_north;
	} else {
		f_copy[slice(d->n * (d->n - 1), d->n, 1)] -= 2 / (d->h_y * d->h_y) * d->boundary_north;
	}
	if (!d->hasNbr(Side::east) && d->neumann) {
		f_copy[slice((d->n - 1), d->n, d->n)] -= 1 / d->h_x * d->boundary_east;
	} else {
		f_copy[slice((d->n - 1), d->n, d->n)] -= 2 / (d->h_x * d->h_x) * d->boundary_east;
	}
	if (!d->hasNbr(Side::south) && d->neumann) {
		f_copy[slice(0, d->n, 1)] += 1 / d->h_y * d->boundary_south;
	} else {
		f_copy[slice(0, d->n, 1)] -= 2 / (d->h_y * d->h_y) * d->boundary_south;
	}
	if (!d->hasNbr(Side::west) && d->neumann) {
		f_copy[slice(0, d->n, d->n)] += 1 / d->h_x * d->boundary_west;
	} else {
		f_copy[slice(0, d->n, d->n)] -= 2 / (d->h_x * d->h_x) * d->boundary_west;
	}

	fftw_execute(plan1);

	tmp /= denom;

	if (d->neumann
	    && !(d->hasNbr(Side::north) || d->hasNbr(Side::east) || d->hasNbr(Side::south)
	         || d->hasNbr(Side::west))) {
		tmp[0] = 0;
	}

	fftw_execute(plan2);

	d->u = u / (4.0 * d->n * d->n);

	if (d->zero_patch) {
		d->u -= d->u.sum() / d->u.size();
	}
}
