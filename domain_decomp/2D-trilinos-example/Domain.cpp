#include "Domain.h"
#include "FftwSolver.h"
#include "FishpackSolver.h"
#include <valarray>
using namespace std;
SolverType Domain::solver_type = SolverType::fftw;
Domain::Domain(DomainSignature ds, int n)
{
	this->n   = n;
	this->h_x = ds.x_length / n;
	this->h_y = ds.y_length / n;
	this->ds  = ds;

#if DD_DEBUG
	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "h:           " << this->h_x << endl;
	cerr << "I start at:  " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "Length:     " << ds.x_length << ", " << ds.y_length << "\n";
	cerr << "North: " << ds.nbr(Side::north) << ", " << ds.nbrRight(Side::north) << "\n";
	cerr << "Idx:   " << ds.index(Side::north) << ", " << ds.indexCenter(Side::north) << "\n";
	cerr << "       " << indexRefinedLeft(Side::north) << ", " << indexRefinedRight(Side::north)
	     << "\n";
	cerr << "East:  " << ds.nbr(Side::east) << ", " << ds.nbrRight(Side::east) << "\n";
	cerr << "Idx:   " << ds.index(Side::east) << ", " << ds.indexCenter(Side::east) << "\n";
	cerr << "       " << indexRefinedLeft(Side::east) << ", " << indexRefinedRight(Side::east)
	     << "\n";
	cerr << "South: " << ds.nbr(Side::south) << ", " << ds.nbrRight(Side::south) << "\n";
	cerr << "Idx:   " << ds.index(Side::south) << ", " << ds.indexCenter(Side::south) << "\n";
	cerr << "       " << indexRefinedLeft(Side::south) << ", " << indexRefinedRight(Side::south)
	     << "\n";
	cerr << "West:  " << ds.nbr(Side::west) << ", " << ds.nbrRight(Side::west) << "\n";
	cerr << "Idx:   " << ds.index(Side::west) << ", " << ds.indexCenter(Side::west) << "\n";
	cerr << "       " << indexRefinedLeft(Side::west) << ", " << indexRefinedRight(Side::west)
	     << "\n";
	cerr << "\n";
#endif
	f      = valarray<double>(n * n);
	f_back = valarray<double>(n * n);
	f_comp = valarray<double>(n * n);
	exact  = valarray<double>(n * n);
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
	delete solver;
#if DD_DEBUG
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

void Domain::plan()
{
	switch (solver_type) {
		case SolverType::fftw:
			solver = new FftwSolver(this);
			break;
		case SolverType::fishpack:
			solver = new FishpackSolver(this);
			break;
	}
}
void Domain::planDirichlet() { plan(); }
void Domain::planNeumann()
{
	neumann = true;
	plan();
}

void Domain::solve() { solver->solve(); }
void Domain::getFluxDiff(vector_type &flux)
{
	auto ptr = flux.getVectorNonConst(0);
	if (hasNbr(Side::north)) {
		fillFluxVector(Side::north, *ptr);
	}
	if (hasNbr(Side::east)) {
		fillFluxVector(Side::east, *ptr);
	}
	if (hasNbr(Side::south)) {
		fillFluxVector(Side::south, *ptr);
	}
	if (hasNbr(Side::west)) {
		fillFluxVector(Side::west, *ptr);
	}
}
void Domain::putCoords(vector_type &xy)
{
	if (hasNbr(Side::north)) {
		if (hasFineNbr(Side::north) || hasCoarseNbr(Side::north)) {
			fillCoords(Side::north, xy);
		}
	}
	if (hasNbr(Side::east)) {
		if (hasFineNbr(Side::east) || hasCoarseNbr(Side::east)) {
			fillCoords(Side::east, xy);
		}
	}
	if (hasNbr(Side::south)) {
		fillCoords(Side::south, xy);
	}
	if (hasNbr(Side::west)) {
		fillCoords(Side::west, xy);
	}
}
void Domain::fillCoords(Side s, vector_type &xy)
{
	auto xy_view = xy.getLocalView<Kokkos::HostSpace>();
	/*if (hasFineNbr(s)) {
	    valarray<double> side        = getSide(s);
	    valarray<double> coarse      = getSideCoarse(s);
	    int              curr_i      = index(s) * n;
	    int              curr_low_i  = indexRefinedRight(s) * n;
	    int              curr_high_i = indexRefinedLeft(s) * n;
	    if (s == Side::north || s == Side::west) {
	        int tmp     = curr_high_i;
	        curr_high_i = curr_low_i;
	        curr_low_i  = tmp;
	    }
	    for (int i = 0; i < n; i++) {
	        diff_view(curr_i, 0) += side[i];
	        diff_view(curr_i, 0) -= 8.0 / 15.0 * (coarse[2 * i] + coarse[2 * i + 1]);
	        diff_view(curr_low_i, 0) += coarse[i] * 8.0 / 15.0;
	        diff_view(curr_high_i, 0) += coarse[n + i] * 8.0 / 15.0;
	        curr_i++;
	        curr_low_i++;
	        curr_high_i++;
	    }
	} else if (hasCoarseNbr(s)) {
	    valarray<double> center        = getSideFine(s);
	    valarray<double> incenter      = getInnerSideFine(s);
	    valarray<double> side          = getSide(s);
	    valarray<double> inside        = getInnerSide(s);
	    int              curr_i        = index(s) * n;
	    int              curr_i_center = indexCenter(s) * n;
	    for (int i = 0; i < n; i++) {
	        diff_view(curr_i_center, 0) += center[i];
	        diff_view(curr_i_center, 0) -= 2.0 / 3.0 * center[i] - 1.0 / 5.0 * incenter[i];
	        diff_view(curr_i, 0) += 2.0 / 3.0 * side[i] - 1.0 / 5.0 * inside[i];
	        curr_i_center++;
	        curr_i++;
	    }
	} else {
	*/
	int curr_i = index(s) * n;
	if (s == Side::north || s == Side::south) {
		double y = ds.y_start;
		if (s == Side::north) {
			y += ds.y_length;
		}
		for (int i = 0; i < n; i++) {
			double x = x_start + h_x / 2.0 + x_length * i / n;
			xy_view(curr_i + i, 0) = x;
			xy_view(curr_i + i, 1) = y;
		}
	} else {
		double x = ds.x_start;
		if (s == Side::east) {
			x += ds.x_length;
		}
		for (int i = 0; i < n; i++) {
			double y = y_start + h_y / 2.0 + y_length * i / n;
			xy_view(curr_i + i, 0) = x;
			xy_view(curr_i + i, 1) = y;
		}
	}
	//}
}
void Domain::fillFluxVector(Side s, single_vector_type &diff)
{
	auto diff_view = diff.getLocalView<Kokkos::HostSpace>();
	if (hasFineNbr(s)) {
		valarray<double> diff = getDiff(s);

		int curr_i = index(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += diff[i];
			curr_i++;
		}
	} else if (hasCoarseNbr(s)) {
		valarray<double> diff_coarse = getDiffCombined(s);

		int curr_i_center = indexCenter(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i_center, 0) += diff_coarse[i];
			curr_i_center++;
		}
	} else {
		valarray<double> diff = getDiff(s);

		int curr_i = index(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += diff[i];
			curr_i++;
		}
	}
}
void Domain::putGhostCells(vector_type &ghost)
{
	auto left_ptr  = ghost.getVectorNonConst(0);
	auto right_ptr = ghost.getVectorNonConst(1);
	if (hasNbr(Side::north)) {
		if (hasFineNbr(Side::north) || hasCoarseNbr(Side::north)) {
			fillGhostVector(Side::north, *left_ptr);
		} else {
			fillGhostVector(Side::north, *right_ptr);
		}
	}
	if (hasNbr(Side::east)) {
		if (hasFineNbr(Side::east) || hasCoarseNbr(Side::east)) {
			fillGhostVector(Side::east, *left_ptr);
		} else {
			fillGhostVector(Side::east, *right_ptr);
		}
	}
	if (hasNbr(Side::south)) {
		fillGhostVector(Side::south, *left_ptr);
	}
	if (hasNbr(Side::west)) {
		fillGhostVector(Side::west, *left_ptr);
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
double Domain::residual()
{
	// integrate in x secton
	double center, north, east, south, west;
	// west
	for (int j = 0; j < n; j++) {
		west   = boundary_west[j];
		center = u[j * n];
		east   = u[j * n + 1];
		if (neumann && !hasNbr(Side::west)) {
			f_comp[j * n] = (-h_x * west - center + east) / (h_x * h_x);
		} else {
			f_comp[j * n] = (2 * west - 3 * center + east) / (h_x * h_x);
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
		} else {
			f_comp[j * n + n - 1] = (west - 3 * center + 2 * east) / (h_x * h_x);
		}
	}
	// south
	for (int i = 0; i < n; i++) {
		south  = boundary_south[i];
		center = u[i];
		north  = u[n + i];
		if (neumann && !hasNbr(Side::south)) {
			f_comp[i] += (-h_y * south - center + north) / (h_y * h_y);
		} else {
			f_comp[i] += (2 * south - 3 * center + north) / (h_y * h_y);
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
		} else {
			f_comp[(n - 1) * n + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
		}
	}
	resid = f - f_comp;
	return sqrt(pow(f - f_comp, 2).sum());
}
void Domain::solveWithInterface(const vector_type &gamma)
{
	auto vec_ptr = gamma.getVector(0);

	Side s = Side::north;
	do {
		if (hasNbr(s)) {
			fillBoundary(s, *vec_ptr);
		}
		s++;
	} while (s != Side::north);

	// solve
	solve();
}
void Domain::getDiff(vector_type &diff)
{
	auto diff_ptr = diff.getVectorNonConst(0);
	Side s        = Side::north;
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
double Domain::integrateU() { return u.sum() * h_x * h_y; }
double Domain::exactNorm() { return sqrt((exact * exact).sum()); }
double Domain::fNorm() { return sqrt((f * f).sum()); }
double Domain::exactNorm(double eavg) { return sqrt(pow(exact - eavg, 2).sum()); }
double                          Domain::integrateExact() { return exact.sum() * h_x * h_y; }
valarray<double> Domain::getSide(const Side s) const
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
valarray<double> Domain::getInnerSide(const Side s) const
{
	valarray<double> retval(n);
	switch (s) {
		case Side::north:
			retval = u[slice(n * (n - 2), n, 1)];
			break;
		case Side::east:
			retval = u[slice(n - 2, n, n)];
			break;
		case Side::south:
			retval = u[slice(n, n, 1)];
			break;
		case Side::west:
			retval = u[slice(1, n, n)];
	}
	return retval;
}
valarray<double> Domain::getSideFineLeft(const Side s) const
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
valarray<double> Domain::getSideFineRight(const Side s) const
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
valarray<double> Domain::getSideFine(const Side s) const
{
	valarray<double> retval(n);
	if (isCoarseLeft(s)) {
		retval = getSideFineLeft(s);
	} else {
		retval = getSideFineRight(s);
	}
	return retval;
}
valarray<double> Domain::getInnerSideFineLeft(const Side s) const
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
valarray<double> Domain::getInnerSideFineRight(const Side s) const
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
valarray<double> Domain::getInnerSideFine(const Side s) const
{
	valarray<double> retval(n);
	if (isCoarseLeft(s)) {
		retval = getInnerSideFineLeft(s);
	} else {
		retval = getInnerSideFineRight(s);
	}
	return retval;
}
valarray<double> Domain::getSideCoarse(const Side s) const
{
	valarray<double> retval = getSide(s);
	valarray<double> refined(2 * n);

	int  grid_i = 1;
	bool right  = false;

	refined[0] = (45.0 * retval[0] - 18.0 * retval[1] + 5.0 * retval[2]) / 32.0;
	refined[1] = (21.0 * retval[0] + 14.0 * retval[1] - 3.0 * retval[2]) / 32.0;
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
	refined[2 * n - 2] = (21.0 * retval[n - 1] + 14.0 * retval[n - 2] - 3.0 * retval[n - 3]) / 32.0;
	refined[2 * n - 1] = (45.0 * retval[n - 1] - 18.0 * retval[n - 2] + 5.0 * retval[n - 3]) / 32.0;
	return refined;
}
valarray<double> Domain::getSideCoarseCombined(const Side s) const
{
	valarray<double> retval = getSide(s);
	valarray<double> refined(n);

	refined[0] = (45.0 * retval[0] - 18.0 * retval[1] + 5.0 * retval[2]) / 32.0;
	refined[0] += (21.0 * retval[0] + 14.0 * retval[1] - 3.0 * retval[2]) / 32.0;
	for (int i = 1; i < n - 1; i++) {
		refined[i] = (5.0 * retval[i - 1] + 30.0 * retval[i] - 3.0 * retval[i + 1]) / 32.0;
		refined[i] += (5.0 * retval[i + 1] + 30.0 * retval[i] - 3.0 * retval[i - 1]) / 32.0;
	}
	refined[n - 1] = (21.0 * retval[n - 1] + 14.0 * retval[n - 2] - 3.0 * retval[n - 3]) / 32.0;
	refined[n - 1] += (45.0 * retval[n - 1] - 18.0 * retval[n - 2] + 5.0 * retval[n - 3]) / 32.0;
	return refined;
}
valarray<double> Domain::getSideCoarseRelativeLeft(const Side s) const
{
	valarray<double> side = getSide(s);
	valarray<double> refined(n);

	int  grid_i = 1;
	bool right  = false;

	refined[0] = (45.0 * side[0] - 18.0 * side[1] + 5.0 * side[2]) / 32.0;
	refined[1] = (21.0 * side[0] + 14.0 * side[1] - 3.0 * side[2]) / 32.0;
	for (int i = 2; i < n; i++) {
		if (right) {
			refined[i]
			= (5.0 * side[grid_i + 1] + 30.0 * side[grid_i] - 3.0 * side[grid_i - 1]) / 32.0;
			right = false;
			grid_i++;
		} else {
			refined[i]
			= (5.0 * side[grid_i - 1] + 30.0 * side[grid_i] - 3.0 * side[grid_i + 1]) / 32.0;
			right = true;
		}
	}
	return refined;
}
valarray<double> Domain::getSideCoarseRelativeRight(const Side s) const
{
	valarray<double> side = getSide(s);
	valarray<double> refined(n);

	int  grid_i = n - 2;
	bool left   = false;

	refined[n - 1] = (45.0 * side[n - 1] - 18.0 * side[n - 2] + 5.0 * side[n - 3]) / 32.0;
	refined[n - 2] = (21.0 * side[n - 1] + 14.0 * side[n - 2] - 3.0 * side[n - 3]) / 32.0;
	for (int i = n - 3; i >= 0; i--) {
		if (left) {
			refined[i]
			= (5.0 * side[grid_i - 1] + 30.0 * side[grid_i] - 3.0 * side[grid_i + 1]) / 32.0;
			left = false;
			grid_i--;
		} else {
			refined[i]
			= (5.0 * side[grid_i + 1] + 30.0 * side[grid_i] - 3.0 * side[grid_i - 1]) / 32.0;
			left = true;
		}
	}
	return refined;
}

valarray<double> Domain::getDiff(const Side s) const { return getBoundary(s) - getSide(s); }
valarray<double> Domain::getDiffFine(const Side s) const
{
	return 2.0 * getBoundary(s)
	       - (getSide(s) + 2.0 / 3.0 * getSide(s) - 1.0 / 5.0 * getInnerSide(s));
}

valarray<double> Domain::getDiffFineToCoarseLeft(const Side s) const
{
	return -1.0 * (getSideFineLeft(s)
	               - (10.0 * getSideFineLeft(s) - 3.0 * getInnerSideFineLeft(s)) / 15.0);
}

valarray<double> Domain::getDiffFineToCoarseRight(const Side s) const
{
	return -1.0 * (getSideFineRight(s)
	               - (10.0 * getSideFineRight(s) - 3.0 * getInnerSideFineRight(s)) / 15.0);
}

valarray<double> Domain::getDiffFineToCoarse(const Side s) const
{
	valarray<double> retval(n);
	if (isCoarseLeft(s)) {
		retval = getDiffFineToCoarseLeft(s);
	} else {
		retval = getDiffFineToCoarseRight(s);
	}
	return retval;
}
valarray<double> Domain::getDiffCombinedLeft(const Side s) const
{
	valarray<double> retval(n);
	valarray<double> side     = getSide(s);
	valarray<double> boundary = getBoundary(s);
	if (s == Side::north || s == Side::west) {
		for (int i = 0; i < n; i++) {
			retval[n - 1 - i / 2] -= side[n - 1 - i];
			retval[n - 1 - i / 2] += boundary[n - 1 - i];
		}
	} else {
		for (int i = 0; i < n; i++) {
			retval[i / 2] -= side[i];
			retval[i / 2] += boundary[i];
		}
	}
	return retval;
}

valarray<double> Domain::getDiffCombinedRight(const Side s) const
{
	valarray<double> retval(n);
	valarray<double> side     = getSide(s);
	valarray<double> boundary = getBoundary(s);
	if (s == Side::north || s == Side::west) {
		for (int i = 0; i < n; i++) {
			retval[i / 2] -= side[i];
			retval[i / 2] += boundary[i];
		}
	} else {
		for (int i = 0; i < n; i++) {
			retval[n - 1 - i / 2] -= side[n - 1 - i];
			retval[n - 1 - i / 2] += boundary[n - 1 - i];
		}
	}
	return retval;
}

valarray<double> Domain::getDiffCombined(const Side s) const
{
	valarray<double> retval(n);
	if (isCoarseLeft(s)) {
		retval = getDiffCombinedLeft(s);
	} else {
		retval = getDiffCombinedRight(s);
	}
	return retval;
}

valarray<double> Domain::getDiffCoarse(const Side s) const
{
	return 2.0 * getBoundary(s)
	       - (getSide(s) + (15.0 * getSide(s) - 8.0 * getSideCoarseCombined(s)) / 15.0);
}
valarray<double> Domain::getDiffCoarseToFineLeft(const Side s) const
{
	return -8.0 / 15.0 * getSideCoarseLeft(s);
}
valarray<double> Domain::getDiffCoarseToFineRight(const Side s) const
{
	return -8.0 / 15.0 * getSideCoarseRight(s);
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

void Domain::fillDiffVector(Side s, single_vector_type &diff)
{
	auto diff_view = diff.getLocalView<Kokkos::HostSpace>();
	if (hasFineNbr(s)) {
		valarray<double> diff       = getDiffCoarse(s);
		valarray<double> diff_left  = getDiffCoarseToFineLeft(s);
		valarray<double> diff_right = getDiffCoarseToFineRight(s);

		int curr_i       = index(s) * n;
		int curr_left_i  = indexRefinedLeft(s) * n;
		int curr_right_i = indexRefinedRight(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += diff[i];
			diff_view(curr_left_i, 0) += diff_left[i];
			diff_view(curr_right_i, 0) += diff_right[i];
			curr_i++;
			curr_left_i++;
			curr_right_i++;
		}
	} else if (hasCoarseNbr(s)) {
		valarray<double> diff        = getDiffFine(s);
		valarray<double> diff_coarse = getDiffFineToCoarse(s);

		int curr_i        = index(s) * n;
		int curr_i_center = indexCenter(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += diff[i];
			diff_view(curr_i_center, 0) += diff_coarse[i];
			curr_i_center++;
			curr_i++;
		}
	} else {
		valarray<double> diff = getDiff(s);

		int curr_i = index(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += diff[i];
			curr_i++;
		}
	}
}
void Domain::fillGhostVector(Side s, single_vector_type &diff)
{
	// weight=false;
	auto diff_view = diff.getLocalView<Kokkos::HostSpace>();
	if (hasFineNbr(s)) {
		valarray<double> side        = getSide(s);
		valarray<double> coarse      = getSideCoarse(s);
		int              curr_i      = index(s) * n;
		int              curr_low_i  = indexRefinedRight(s) * n;
		int              curr_high_i = indexRefinedLeft(s) * n;
		if (s == Side::north || s == Side::west) {
			int tmp     = curr_high_i;
			curr_high_i = curr_low_i;
			curr_low_i  = tmp;
		}
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) += side[i];
			diff_view(curr_i, 0) -= 8.0 / 15.0 * (coarse[2 * i] + coarse[2 * i + 1]);
			diff_view(curr_low_i, 0) += coarse[i] * 8.0 / 15.0;
			diff_view(curr_high_i, 0) += coarse[n + i] * 8.0 / 15.0;
			curr_i++;
			curr_low_i++;
			curr_high_i++;
		}
	} else if (hasCoarseNbr(s)) {
		valarray<double> center        = getSideFine(s);
		valarray<double> incenter      = getInnerSideFine(s);
		valarray<double> side          = getSide(s);
		valarray<double> inside        = getInnerSide(s);
		int              curr_i        = index(s) * n;
		int              curr_i_center = indexCenter(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i_center, 0) += center[i];
			diff_view(curr_i_center, 0) -= 2.0 / 3.0 * center[i] - 1.0 / 5.0 * incenter[i];
			diff_view(curr_i, 0) += 2.0 / 3.0 * side[i] - 1.0 / 5.0 * inside[i];
			curr_i_center++;
			curr_i++;
		}
	} else {
		valarray<double> side   = getSide(s);
		int              curr_i = index(s) * n;
		for (int i = 0; i < n; i++) {
			diff_view(curr_i, 0) = side[i];
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
double Domain::area() { return x_length * y_length; }
double Domain::integrateBoundaryFlux()
{
	double sum = 0;
	if (!hasNbr(Side::north)) {
		sum += boundary_north.sum() * h_x;
	}
	if (!hasNbr(Side::east)) {
		sum -= boundary_east.sum() * h_y;
	}
	if (!hasNbr(Side::south)) {
		sum -= boundary_south.sum() * h_x;
	}
	if (!hasNbr(Side::west)) {
		sum += boundary_west.sum() * h_y;
	}
	return sum;
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
			os << u[loc] << tab << error[loc] << tab << resid[loc] * h_x * h_y << endl;
		}
		os << endl;
	}
}
