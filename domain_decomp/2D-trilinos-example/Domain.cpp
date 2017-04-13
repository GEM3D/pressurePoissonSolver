#include "Domain.h"
#include "Iface.h"
using namespace std;
Domain::Domain(DomainSignature ds, int n, double h_x, double h_y)
{
	this->n  = n;
	this->h_x = h_x / ds.refine_level;
	this->h_y = h_y / ds.refine_level;
	this->ds = ds;

#if NDEBUG
	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "h:           " << this->h_x << endl;
	cerr << "I start at:  " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "Length:     " << ds.x_length << ", " << ds.y_length << "\n";
	cerr << "North: " << ds.nbr(Side::north) << ", " << ds.nbrRight(Side::north) << "\n";
	cerr << "Idx:   " << globalIndex(Side::north) << ", " << globalIndexCenter(Side::north) << "\n";
	cerr << "East:  " << ds.nbr(Side::east) << ", " << ds.nbrRight(Side::east) << "\n";
	cerr << "Idx:   " << globalIndex(Side::east) << ", " << globalIndexCenter(Side::east) << "\n";
	cerr << "South: " << ds.nbr(Side::south) << ", " << ds.nbrRight(Side::south) << "\n";
	cerr << "Idx:   " << globalIndex(Side::south) << ", " << globalIndexCenter(Side::south) << "\n";
	cerr << "West:  " << ds.nbr(Side::west) << ", " << ds.nbrRight(Side::west) << "\n";
	cerr << "Idx:   " << globalIndex(Side::west) << ", " << globalIndexCenter(Side::west) << "\n";
	cerr << "\n";
#endif
	f      = valarray<double>(n * n);
	f_back = valarray<double>(n * n);
	f_copy = valarray<double>(n * n);
	exact  = valarray<double>(n * n);
	tmp    = valarray<double>(n * n);
	u      = valarray<double>(n * n);
	u_back = valarray<double>(n * n);
	denom  = valarray<double>(n * n);

	boundary_north = valarray<double>(n);
	boundary_south = valarray<double>(n);
	boundary_east  = valarray<double>(n);
	boundary_west  = valarray<double>(n);

	boundary_north_back = valarray<double>(n);
	boundary_south_back = valarray<double>(n);
	boundary_east_back  = valarray<double>(n);
	boundary_west_back  = valarray<double>(n);

	x_start  = ds.x_start;
	y_start  = ds.y_start;
	x_length = ds.x_length;
	y_length = ds.y_length;
}

Domain::~Domain()
{
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
#if NDEBUG
	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "h:           " << this->h_x << endl;
	cerr << "I start at:  " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "Length:     " << ds.x_length << ", " << ds.y_length << "\n";
	cerr << "North: " << ds.nbr(Side::north) << ", " << ds.nbrRight(Side::north) << "\n";
	cerr << "Idx:   " << index(Side::north) << ", " << indexCenter(Side::north) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::north) << ", " << hasCoarseNbr(Side::north) << ", "
	     << isCoarseLeft(Side::north) << endl;
	cerr << "East:  " << ds.nbr(Side::east) << ", " << ds.nbrRight(Side::east) << "\n";
	cerr << "Idx:   " << index(Side::east) << ", " << indexCenter(Side::east) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::east) << ", " << hasCoarseNbr(Side::east) << ", "
	     << isCoarseLeft(Side::east) << endl;
	cerr << "South: " << ds.nbr(Side::south) << ", " << ds.nbrRight(Side::south) << "\n";
	cerr << "Idx:   " << index(Side::south) << ", " << indexCenter(Side::south) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::south) << ", " << hasCoarseNbr(Side::south) << ", "
	     << isCoarseLeft(Side::south) << endl;
	cerr << "West:  " << ds.nbr(Side::west) << ", " << ds.nbrRight(Side::west) << "\n";
	cerr << "Idx:   " << index(Side::west) << ", " << indexCenter(Side::west) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::west) << ", " << hasCoarseNbr(Side::west) << ", "
	     << isCoarseLeft(Side::west) << endl;
	cerr << "\n";
#endif
}

void Domain::planDirichlet()
{
	// create plan
	plan1 = fftw_plan_r2r_2d(n, n, &f_copy[0], &tmp[0], FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

	plan2 = fftw_plan_r2r_2d(n, n, &tmp[0], &u[0], FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);
	// create denom vector

	for (int yi = 0; yi < n; yi++) {
		denom[slice(yi * n, n, 1)] = -4 / (h_x * h_x) * pow(sin((yi + 1) * M_PI / (2 * n)), 2);
	}

	valarray<double> ones(n);
	ones = 1;
	for (int xi = 0; xi < n; xi++) {
		denom[slice(xi, n, n)]
		-= 4 / (h_y * h_y) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
	};

}

void Domain::planNeumann()
{
	neumann                       = true;
	fftw_r2r_kind x_transform     = FFTW_RODFT10;
	fftw_r2r_kind x_transform_inv = FFTW_RODFT01;
	fftw_r2r_kind y_transform     = FFTW_RODFT10;
	fftw_r2r_kind y_transform_inv = FFTW_RODFT01;
	if (!hasNbr(Side::east) && !hasNbr(Side::west)) {
		x_transform     = FFTW_REDFT10;
		x_transform_inv = FFTW_REDFT01;
	} else if (!hasNbr(Side::west)) {
		x_transform     = FFTW_REDFT11;
		x_transform_inv = FFTW_REDFT11;
	} else if (!hasNbr(Side::east)) {
		x_transform     = FFTW_RODFT11;
		x_transform_inv = FFTW_RODFT11;
	}
	if (!hasNbr(Side::north) && !hasNbr(Side::south)) {
		y_transform     = FFTW_REDFT10;
		y_transform_inv = FFTW_REDFT01;
	} else if (!hasNbr(Side::south)) {
		y_transform     = FFTW_REDFT11;
		y_transform_inv = FFTW_REDFT11;
	} else if (!hasNbr(Side::north)) {
		y_transform     = FFTW_RODFT11;
		y_transform_inv = FFTW_RODFT11;
	}
	plan1 = fftw_plan_r2r_2d(n, n, &f_copy[0], &tmp[0], y_transform, x_transform, FFTW_MEASURE);

	plan2
	= fftw_plan_r2r_2d(n, n, &tmp[0], &u[0], y_transform_inv, x_transform_inv, FFTW_MEASURE);

	// create denom vector
	if (!hasNbr(Side::north) && !hasNbr(Side::south)) {
		for (int yi = 0; yi < n; yi++) {
			denom[slice(yi * n, n, 1)] = -4 / (h_x * h_x) * pow(sin(yi * M_PI / (2 * n)), 2);
		}
	} else if (!hasNbr(Side::south) || !hasNbr(Side::north)) {
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

	valarray<double> ones(n);
	ones = 1;

	if (!hasNbr(Side::east) && !hasNbr(Side::west)) {
		for (int xi = 0; xi < n; xi++) {
			denom[slice(xi, n, n)] -= 4 / (h_y * h_y) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
		}
	} else if (!hasNbr(Side::west) || !hasNbr(Side::east)) {
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

void Domain::solve()
{
	f_copy = f;
	if (!hasNbr(Side::north) && neumann) {
		f_copy[slice(n * (n - 1), n, 1)] -= 1 / h_y * boundary_north;
	} else {
		f_copy[slice(n * (n - 1), n, 1)] -= 2 / (h_y * h_y) * boundary_north;
	}
	if (!hasNbr(Side::east) && neumann) {
		f_copy[slice((n - 1), n, n)] -= 1 / h_x * boundary_east;
	} else {
		f_copy[slice((n - 1), n, n)] -= 2 / (h_x * h_x) * boundary_east;
	}
	if (!hasNbr(Side::south) && neumann) {
		f_copy[slice(0, n, 1)] += 1 / h_y * boundary_south;
	} else {
		f_copy[slice(0, n, 1)] -= 2 / (h_y * h_y) * boundary_south;
	}
	if (!hasNbr(Side::west) && neumann) {
		f_copy[slice(0, n, n)] += 1 / h_x * boundary_west;
	} else {
		f_copy[slice(0, n, n)] -= 2 / (h_x * h_x) * boundary_west;
	}

	fftw_execute(plan1);

	tmp /= denom;

	if (neumann
	    && !(hasNbr(Side::north) || hasNbr(Side::east) || hasNbr(Side::south)
	         || hasNbr(Side::west))) {
		tmp[0] = 0;
	}

	fftw_execute(plan2);

	u /= 4 * n * n;
}

void Domain::putGhostCells(vector_type &ghost)
{
	auto left_ptr  = ghost.getVectorNonConst(0);
	auto right_ptr = ghost.getVectorNonConst(1);
	if (hasNbr(Side::north)) {
		if (hasFineNbr(Side::north) || hasCoarseNbr(Side::north)) {
			fillDiffVector(Side::north, *left_ptr, true);
		} else {
			fillDiffVector(Side::north, *right_ptr, true);
		}
	}
	if (hasNbr(Side::east)) {
		if (hasFineNbr(Side::east) || hasCoarseNbr(Side::east)) {
			fillDiffVector(Side::east, *left_ptr, true);
		} else {
			fillDiffVector(Side::east, *right_ptr, true);
		}
	}
	if (hasNbr(Side::south)) {
		fillDiffVector(Side::south, *left_ptr, true);
	};
	if (hasNbr(Side::west)) {
		fillDiffVector(Side::west, *left_ptr, true);
	}
}

double Domain::residual(vector_type &ghost)
{
	auto left_ptr  = ghost.getVector(0);
	auto right_ptr = ghost.getVector(1);

	if (hasNbr(Side::north)) {
	    fillBoundary(Side::north, *left_ptr);
	}
	if (hasNbr(Side::east)) {
	    fillBoundary(Side::east, *left_ptr);
	}
	if (hasNbr(Side::south)) {
		if (hasFineNbr(Side::south) || hasCoarseNbr(Side::south)) {
			fillBoundary(Side::south, *left_ptr);
		} else {
			fillBoundary(Side::south, *right_ptr);
		}
	}
	if (hasNbr(Side::west)) {
		if (hasFineNbr(Side::west) || hasCoarseNbr(Side::west)) {
			fillBoundary(Side::west, *left_ptr);
		} else {
			fillBoundary(Side::west, *right_ptr);
		}
	}

	valarray<double> f_comp = valarray<double>(n * n);
	// integrate in x secton
	double center, north, east, south, west;
	// west
	for (int j = 0; j < n; j++) {
		west   = boundary_west[j];
		center = u[j * n];
		east   = u[j * n + 1];
		if (neumann && !hasNbr(Side::west)) {
			f_comp[j * n] = (-h_x * west - center + east) / (h_x * h_x);
		} else if (!hasNbr(Side::west)) {
			f_comp[j * n] = (2 * west - 3 * center + east) / (h_x * h_x);
		} else {
		    f_comp[j * n] = (west - 2 * center + east) / (h_x * h_x);
		}
	}
	// middle
	for (int i = 1; i < n - 1; i++) {
		for (int j = 0; j < n; j++) {
			east   = u[j * n + i - 1];
			center = u[j * n + i];
			west   = u[j * n + i + 1];

			f_comp[j * n + i] = (west - 2 * center + east) / (h_x * h_x);
		}
	}
	// east
	for (int j = 0; j < n; j++) {
		west   = u[j * n + n - 2];
		center = u[j * n + n - 1];
		east   = boundary_east[j];
		if (neumann && !hasNbr(Side::east)) {
			f_comp[j * n + n - 1] = (west - center + h_x * east) / (h_x * h_x);
		} else if (!hasNbr(Side::east)) {
			f_comp[j * n + n - 1] = (west - 3 * center + 2 * east) / (h_x * h_x);
		} else {
			f_comp[j * n + n - 1] = (west - 2 * center + east) / (h_x * h_x);
		}
	}
	// south
	for (int i = 0; i < n; i++) {
		south  = boundary_south[i];
		center = u[i];
		north  = u[n + i];
		if (neumann && !hasNbr(Side::south)) {
			f_comp[i] += (-h_y * south - center + north) / (h_y * h_y);
		} else if (!hasNbr(Side::south)) {
			f_comp[i] += (2 * south - 3 * center + north) / (h_y * h_y);
		} else {
			f_comp[i] += (south - 2 * center + north) / (h_y * h_y);
		}
	}
	// middle
	for (int i = 0; i < n; i++) {
		for (int j = 1; j < n - 1; j++) {
			south  = u[(j - 1) * n + i];
			center = u[j * n + i];
			north  = u[(j + 1) * n + i];

			f_comp[j * n + i] += (south - 2 * center + north) / (h_y * h_y);
		}
	}
	// north
	for (int i = 0; i < n; i++) {
		south  = u[(n - 2) * n + i];
		center = u[(n - 1) * n + i];
		north  = boundary_north[i];
		if (neumann && !hasNbr(Side::north)) {
			f_comp[(n - 1) * n + i] += (south - center + h_y * north) / (h_y * h_y);
		} else if (!hasNbr(Side::north)) {
			f_comp[(n - 1) * n + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
		} else {
			f_comp[(n - 1) * n + i] += (south - 2 * center + north) / (h_y * h_y);
		}
	}
	resid = f - f_comp;
	return sqrt(pow(f - f_comp, 2).sum());
}

void Domain::solveWithInterface(const vector_type &gamma, vector_type &diff)
{
	auto vec_ptr  = gamma.getVector(0);
	auto diff_ptr = diff.getVectorNonConst(0);

	Side s = Side::north;
	do {
		if (hasNbr(s)) {
			fillBoundary(s, *vec_ptr);
		}
		s++;
	} while (s != Side::north);

	// solve
	solve();

	s = Side::north;
	do {
		if (hasNbr(s)) {
			fillDiffVector(s, *diff_ptr);
		}
		s++;
	} while (s != Side::north);
}
double Domain::diffNorm()
{
	error = exact - u;
	return sqrt(pow(exact - u, 2).sum());
}
double Domain::diffNorm(double uavg, double eavg)
{
	error = exact - u - eavg + uavg;
	return sqrt(pow(exact - u - eavg + uavg, 2).sum());
}
double Domain::uSum() { return u.sum(); }
double Domain::exactNorm() { return sqrt((exact * exact).sum()); }
double Domain::fNorm() { return sqrt((f * f).sum()); }
double Domain::exactNorm(double eavg) { return sqrt(pow(exact - eavg, 2).sum()); }
double                          Domain::exactSum() { return exact.sum(); }
valarray<double> Domain::getSide(Side s)
{
	valarray<double> retval(n);
	switch (s) {
		case Side::north:
			retval = u[slice(n * (n - 1), n, 1)];
			break;
		case Side::east:
			retval = u[slice(n - 1, n, n)];
			break;
		case Side::south:
			retval = u[slice(0, n, 1)];
			break;
		case Side::west:
			retval = u[slice(0, n, n)];
	}
	return retval;
}
valarray<double> Domain::getInnerSide(Side s)
{
	valarray<double> retval(n);
	switch (s) {
		case Side::north:
			retval = u[slice(n * (n - 1), n, 1)];
			break;
		case Side::east:
			retval = u[slice(n - 1, n, n)];
			break;
		case Side::south:
			retval = u[slice(0, n, 1)];
			break;
		case Side::west:
			retval = u[slice(1, n, n)];
	}
	return retval;
}
valarray<double> Domain::getSideFineLeft(Side s)
{
	valarray<double> retval(n);
	valarray<double> side = getSide(s);
	if (s == Side::north || s == Side::west) {
		for (int i = 0; i < n; i++) {
			retval[n - 1 - i / 2] += side[n - 1 - i];
		}
	} else {
		for (int i = 0; i < n; i++) {
			retval[i / 2] += side[i];
		}
	}
	return retval;
}
valarray<double> Domain::getSideFineRight(Side s)
{
	valarray<double> retval(n);
	valarray<double> side = getSide(s);
	if (s == Side::north || s == Side::west) {
		for (int i = 0; i < n; i++) {
			retval[i / 2] += side[i];
		}
	} else {
		for (int i = 0; i < n; i++) {
			retval[n - 1 - i / 2] += side[n - 1 - i];
		}
	}
	return retval;

}
valarray<double> Domain::getSideFine(Side s)
{
	valarray<double> retval(n);
	if (isCoarseLeft(s)) {
        retval = getSideFineLeft(s);
	} else {
        retval = getSideFineRight(s);
	}
	return retval;
}
valarray<double> Domain::getInnerSideFineLeft(Side s)
{
	valarray<double> retval(n);
	valarray<double> side = getInnerSide(s);
	if (s == Side::north || s == Side::west) {
		for (int i = 0; i < n; i++) {
			retval[n - 1 - i / 2] += side[n - 1 - i];
		}
	} else {
		for (int i = 0; i < n; i++) {
			retval[i / 2] += side[i];
		}
	}
	return retval;
}
valarray<double> Domain::getInnerSideFineRight(Side s)
{
	valarray<double> retval(n);
	valarray<double> side = getInnerSide(s);
	if (s == Side::north || s == Side::west) {
		for (int i = 0; i < n; i++) {
			retval[i / 2] += side[i];
		}
	} else {
		for (int i = 0; i < n; i++) {
			retval[n - 1 - i / 2] += side[n - 1 - i];
		}
	}
	return retval;

}
valarray<double> Domain::getInnerSideFine(Side s)
{
	valarray<double> retval(n);
	if (isCoarseLeft(s)) {
        retval = getInnerSideFineLeft(s);
	} else {
        retval = getInnerSideFineRight(s);
	}
	return retval;
}
valarray<double> Domain::getSideCoarse(Side s)
{
	bool   nbr_left=false;
	bool   nbr_right=false;
	double val_left=0;
	double val_right=0;
	switch (s) {
		case Side::north:
			nbr_left  = hasNbr(Side::west);
			val_left  = boundary_west[n - 1];
			nbr_right = hasNbr(Side::east);
			val_right = boundary_east[n - 1];
			break;
		case Side::east:
			nbr_left  = hasNbr(Side::south);
			val_left  = boundary_south[n - 1];
			nbr_right = hasNbr(Side::north);
			val_right = boundary_north[n - 1];
			break;
		case Side::south:
			nbr_left  = hasNbr(Side::west);
			val_left  = boundary_west[0];
			nbr_right = hasNbr(Side::east);
			val_right = boundary_east[0];
			break;
		case Side::west:
			nbr_left  = hasNbr(Side::south);
			val_left  = boundary_south[0];
			nbr_right = hasNbr(Side::north);
			val_right = boundary_north[0];
	}
    valarray<double> retval = getSide(s);
	valarray<double> refined(2 * n);

	int  grid_i = 1;
	bool right  = false;

	if (!nbr_left && neumann) {
		refined[0] = (-10.0 * val_left + 35.0 * retval[0] - 3.0 * retval[1]) / 32.0;
		refined[1] = (6.0 * val_left + 27.0 * retval[0] + 5.0 * retval[1]) / 32.0;
	} else if (!nbr_left) {
		refined[0] = (10.0 * val_left + 15.0 * retval[0] - retval[1]) / 24.0;
		refined[1] = (-2.0 * val_left + 9.0 * retval[0] + retval[1]) / 8.0;
    }else{
		refined[0] = (10.0 * val_left + 25.0 * retval[0] -3.0* retval[1]) / 32.0;
		refined[1] = (-6.0 * val_left + 33.0 * retval[0] +5.0* retval[1]) / 32.0;
    }
	for (int i = 2; i < 2 * n - 2; i++) {
		if (right) {
			refined[i]
			= (5.0 * retval[grid_i + 1] + 30.0 * retval[grid_i] - 3.0 * retval[grid_i - 1]) / 32.0;
			right = false;
			grid_i++;
		} else {
			refined[i]
			= (5.0 * retval[grid_i - 1] + 30.0 * retval[grid_i] - 3.0 * retval[grid_i + 1]) / 32.0;
			right = true;
		}
	}
	if (!nbr_right && neumann) {
		refined[2 * n - 2] = (6.0 * val_right + 27.0 * retval[n - 1] + 5.0 * retval[n - 2]) / 32.0;
		refined[2 * n - 1]
		= (-10.0 * val_right + 35.0 * retval[n - 1] - 3.0 * retval[n - 2]) / 32.0;
	} else if (!nbr_right) {
		refined[2 * n - 2] = (-2.0 * val_right + 9.0 * retval[n - 1] + retval[n - 2]) / 8.0;
		refined[2 * n - 1] = (10.0 * val_right + 15.0 * retval[n - 1] - retval[n - 2]) / 24.0;
	} else {
		refined[2 * n - 2] = (-6.0 * val_right + 33.0 * retval[n - 1] + 5.0 * retval[n - 2]) / 32.0;
		refined[2 * n - 1] = (10.0 * val_right + 25.0 * retval[n - 1] - 3.0 * retval[n - 2]) / 32.0;
	}
	return refined;
}
valarray<double> Domain::getStencil(Side s,Tilt t)
{
		return getSide(s);
}

void Domain::fillBoundary(Side s, const single_vector_type &gamma)
{
	int               curr_i     = index(s) * n;
	valarray<double> &boundary   = *getBoundaryPtr(s);
	auto              gamma_view = gamma.getLocalView<Kokkos::HostSpace>();

	for (int i = 0; i < n; i++) {
		boundary[i] = gamma_view(curr_i, 0);
		curr_i++;
	}
}
void Domain::fillDiffVector(Side s, single_vector_type &diff, bool residual)
{
	// weight=false;
	auto diff_view = diff.getLocalView<Kokkos::HostSpace>();
	if (hasFineNbr(s)) {
		valarray<double> side   = getSide(s);
		valarray<double> coarse = getSideCoarse(s);
		if (!residual) {
			int curr_i = index(s) * n;
			for (int i = 0; i < n; i++) {
				diff_view(curr_i, 0) += side[i];
				curr_i++;
			}
		}
		int curr_i      = index(s) * n;
		int curr_low_i  = indexRefinedRight(s) * n;
		int curr_high_i = indexRefinedLeft(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += side[i];
			diff_view(curr_i, 0) -= 8.0 / 15.0 * (coarse[2 * i] + coarse[2 * i + 1]);
			diff_view(curr_low_i, 0) += coarse[i] * 8.0 / 15.0;
			diff_view(curr_high_i, 0) += coarse[n + i]*8.0/15.0;
			curr_i++;
			curr_low_i++;
			curr_high_i++;
		}
	} else if (hasCoarseNbr(s)) {
		valarray<double> center = getSideFine(s);
		valarray<double> incenter = getInnerSideFine(s);
		valarray<double> side   = getSide(s);
        valarray<double> inside = getInnerSide(s);
		if (!residual) {
			int curr_i        = index(s) * n;
			int curr_i_center = indexCenter(s) * n;
			for (int i = 0; i < n; i++) {
				diff_view(curr_i_center, 0) += center[i];
				diff_view(curr_i, 0) += side[i];
				curr_i_center++;
				curr_i++;
			}
		}
		int curr_i        = index(s) * n;
		int curr_i_center = indexCenter(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i_center, 0) -= 2.0 / 3.0 * center[i] - 1.0 / 5.0 * incenter[i];
			diff_view(curr_i, 0) += 2.0 / 3.0 * side[i] - 1.0 / 5.0 * inside[i];
			curr_i_center++;
			curr_i++;
		}
	} else {
		valarray<double> side   = getSide(s);
		int              curr_i = index(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += side[i];
			curr_i++;
		}
	}
}
void Domain::swapResidSol()
{
	boundary_north_back = boundary_north;
	boundary_east_back  = boundary_east;
	boundary_south_back = boundary_south;
	boundary_west_back  = boundary_west;
	boundary_north      = 0;
	boundary_east       = 0;
	boundary_south      = 0;
	boundary_west       = 0;
	u_back              = u;
	f_back              = f;
	f                   = resid;
}
void Domain::sumResidIntoSol()
{
	boundary_north = boundary_north_back;
	boundary_east  = boundary_east_back;
	boundary_south = boundary_south_back;
	boundary_west  = boundary_west_back;
	f              = f_back;
	u += u_back;
}
void Domain::outputClaw(std::ostream &os)
{
	const string tab = "\t";
	os << ds.id << tab << "grid_number" << endl;
	os << ds.refine_level << tab << "AMR_level" << endl;
	os << 0 << tab << "block_number" << endl;
	os << 0 << tab << "mpi_rank" << endl;
	os << n << tab << "mx" << endl;
	os << n << tab << "my" << endl;
	os << x_start << tab << "xlow" << endl;
	os << y_start << tab << "ylow" << endl;
	os << h_x << tab << "dx" << endl;
	os << h_y << tab << "dy" << endl;
	os << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int loc = j + i * n;
			os << u[loc] << tab << error[loc] << tab << resid[loc] << endl;
		}
		os << endl;
	}
}
