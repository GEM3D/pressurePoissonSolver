#include "Domain.h"
using namespace std;
Domain::Domain(DomainSignature ds, int nx, int ny, double h_x, double h_y)
{
	this->nx  = nx;
	this->ny  = ny;
	this->h_x = h_x / ds.refine_level;
	this->h_y = h_y / ds.refine_level;

	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "I start at: " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "North: " << ds.nbr[0] << ", " << ds.nbr[1] << "\n";
	cerr << "Idx:   " << ds.global_i[0] << ", " << ds.global_i[1] << "\n";
	cerr << "East:  " << ds.nbr[2] << ", " << ds.nbr[3] << "\n";
	cerr << "Idx:   " << ds.global_i[2] << ", " << ds.global_i[3] << "\n";
	cerr << "South: " << ds.nbr[4] << ", " << ds.nbr[5] << "\n";
	cerr << "Idx:   " << ds.global_i[4] << ", " << ds.global_i[5] << "\n";
	cerr << "West:  " << ds.nbr[6] << ", " << ds.nbr[7] << "\n";
	cerr << "Idx:   " << ds.global_i[6] << ", " << ds.global_i[7] << "\n";
	cerr << "\n";
	f      = valarray<double>(nx * ny);
	f_copy = valarray<double>(nx * ny);
	exact  = valarray<double>(nx * ny);
	tmp    = valarray<double>(nx * ny);
	u      = valarray<double>(nx * ny);
	denom  = valarray<double>(nx * ny);

	boundary_north = valarray<double>(nx);
	boundary_south = valarray<double>(nx);
	boundary_east  = valarray<double>(ny);
	boundary_west  = valarray<double>(ny);

	for (int q = 0; q < 8; q++) {
		if (ds.nbr[q] != -1) {
			global_i[q] = ds.global_i[q] * nx;
			iface_i[q]  = ds.global_i[q] * 22;
		}
	}

	nbr = ds.nbr;

	x_start  = ds.x_start;
	y_start  = ds.y_start;
	x_length = ds.x_length;
	y_length = ds.y_length;
}

Domain::~Domain()
{
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
}

void Domain::planDirichlet()
{
	// create plan
	plan1 = fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

	plan2 = fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);
	// create denom vector

	for (int yi = 0; yi < nx; yi++) {
		denom[slice(yi * ny, nx, 1)]
		= -4 / (h_x * h_x) * pow(sin((yi + 1) * M_PI / (2 * nx)), 2);
	}

	valarray<double> ones(ny);
	ones = 1;
	for (int xi = 0; xi < ny; xi++) {
		denom[slice(xi, ny, nx)]
		-= 4 / (h_y * h_y) * pow(sin((xi + 1) * M_PI / (2 * ny)), 2) * ones;
	}
}

void Domain::planNeumann()
{
	neumann                       = true;
	fftw_r2r_kind x_transform     = FFTW_RODFT10;
	fftw_r2r_kind x_transform_inv = FFTW_RODFT01;
	fftw_r2r_kind y_transform     = FFTW_RODFT10;
	fftw_r2r_kind y_transform_inv = FFTW_RODFT01;
	if (nbr[2] == -1 && nbr[6] == -1) {
		x_transform     = FFTW_REDFT10;
		x_transform_inv = FFTW_REDFT01;
	} else if (nbr[6] == -1) {
		x_transform     = FFTW_REDFT11;
		x_transform_inv = FFTW_REDFT11;
	} else if (nbr[2] == -1) {
		x_transform     = FFTW_RODFT11;
		x_transform_inv = FFTW_RODFT11;
	}
	if (nbr[0] == -1 && nbr[4] == -1) {
		y_transform     = FFTW_REDFT10;
		y_transform_inv = FFTW_REDFT01;
	} else if (nbr[4] == -1) {
		y_transform     = FFTW_REDFT11;
		y_transform_inv = FFTW_REDFT11;
	} else if (nbr[0] == -1) {
		y_transform     = FFTW_RODFT11;
		y_transform_inv = FFTW_RODFT11;
	}
	plan1 = fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], y_transform, x_transform, FFTW_MEASURE);

	plan2
	= fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], y_transform_inv, x_transform_inv, FFTW_MEASURE);

	// create denom vector
	if (nbr[0] == -1 && nbr[4] == -1) {
		for (int yi = 0; yi < nx; yi++) {
			denom[slice(yi * ny, nx, 1)]
			= -4 / (h_x * h_x) * pow(sin(yi * M_PI / (2 * nx)), 2);
		}
	} else if (nbr[4] == -1 || nbr[0] == -1) {
		for (int yi = 0; yi < nx; yi++) {
			denom[slice(yi * ny, nx, 1)]
			= -4 / (h_x * h_x) * pow(sin((yi + 0.5) * M_PI / (2 * nx)), 2);
		}
	} else {
		for (int yi = 0; yi < nx; yi++) {
			denom[slice(yi * ny, nx, 1)]
			= -4 / (h_x * h_x) * pow(sin((yi + 1) * M_PI / (2 * nx)), 2);
		}
	}

	valarray<double> ones(ny);
	ones = 1;

	if (nbr[2] == -1 && nbr[6] == -1) {
		for (int xi = 0; xi < ny; xi++) {
			denom[slice(xi, ny, nx)]
			-= 4 / (h_y * h_y) * pow(sin(xi * M_PI / (2 * ny)), 2) * ones;
		}
	} else if (nbr[6] == -1 || nbr[2] == -1) {
		for (int xi = 0; xi < ny; xi++) {
			denom[slice(xi, ny, nx)]
			-= 4 / (h_y * h_y) * pow(sin((xi + 0.5) * M_PI / (2 * ny)), 2) * ones;
		}
	} else {
		for (int xi = 0; xi < ny; xi++) {
			denom[slice(xi, ny, nx)]
			-= 4 / (h_y * h_y) * pow(sin((xi + 1) * M_PI / (2 * ny)), 2) * ones;
		}
	}
}

void Domain::solve()
{
	f_copy = f;
	if (nbr[0] == -1 && neumann) {
		f_copy[slice(nx * (ny - 1), nx, 1)] -= 1 / h_y * boundary_north;
	} else {
		f_copy[slice(nx * (ny - 1), nx, 1)] -= 2 / (h_y * h_y) * boundary_north;
	}
	if (nbr[2] == -1 && neumann) {
		f_copy[slice((nx - 1), ny, nx)] -= 1 / h_x * boundary_east;
	} else {
		f_copy[slice((nx - 1), ny, nx)] -= 2 / (h_x * h_x) * boundary_east;
	}
	if (nbr[4] == -1 && neumann) {
		f_copy[slice(0, nx, 1)] += 1 / h_y * boundary_south;
	} else {
		f_copy[slice(0, nx, 1)] -= 2 / (h_y * h_y) * boundary_south;
	}
	if (nbr[6] == -1 && neumann) {
		f_copy[slice(0, ny, nx)] += 1 / h_x * boundary_west;
	} else {
		f_copy[slice(0, ny, nx)] -= 2 / (h_x * h_x) * boundary_west;
	}

	fftw_execute(plan1);

	tmp /= denom;

	if (neumann && nbr[0] == -1 && nbr[2] == -1 && nbr[4] == -1 && nbr[6] == -1) {
		tmp[0] = 0;
	}

	fftw_execute(plan2);

	u /= 4 * nx * ny;
}

void Domain::putGhostCells(vector_type &ghost)
{
	auto ghost_view = ghost.getLocalView<Kokkos::HostSpace>();
	if (nbr[0] != -1) {
		int curr_i = local_i[0];
		for (int i = 0; i < nx; i++) {
			ghost_view(curr_i, 1) = u[nx * (ny - 1) + i];
			curr_i++;
		}
	}
	if (nbr[2] != -1) {
		int curr_i = local_i[2];
		for (int i = 0; i < ny; i++) {
			ghost_view(curr_i, 1) = u[(i + 1) * nx - 1];
			curr_i++;
		}
	}
	if (nbr[4] != -1) {
		int curr_i = local_i[4];
		for (int i = 0; i < nx; i++) {
			ghost_view(curr_i, 0) = u[i];
			curr_i++;
		}
	}
	if (nbr[6] != -1) {
		int curr_i = local_i[6];
		for (int i = 0; i < ny; i++) {
			ghost_view(curr_i, 0) = u[i * nx];
			curr_i++;
		}
	}
}

double Domain::residual(vector_type &ghost)
{
	auto ghost_view = ghost.getLocalView<Kokkos::HostSpace>();
	if (nbr[0] != -1) {
		int curr_i = local_i[0];
		for (int i = 0; i < nx; i++) {
			boundary_north[i] = ghost_view(curr_i, 0);
			curr_i++;
		}
	}
	if (nbr[2] != -1) {
		int curr_i = local_i[2];
		for (int i = 0; i < ny; i++) {
			boundary_east[i] = ghost_view(curr_i, 0);
			curr_i++;
		}
	}
	if (nbr[4] != -1) {
		int curr_i = local_i[4];
		for (int i = 0; i < nx; i++) {
			boundary_south[i] = ghost_view(curr_i, 1);
			curr_i++;
		}
	}
	if (nbr[6] != -1) {
		int curr_i = local_i[6];
		for (int i = 0; i < ny; i++) {
			boundary_west[i] = ghost_view(curr_i, 1);
			curr_i++;
		}
	}

	valarray<double> f_comp = valarray<double>(nx * ny);
	// integrate in x directon
	double center, north, east, south, west;
	// west
	for (int j = 0; j < ny; j++) {
		west   = boundary_west[j];
		center = u[j * nx];
		east   = u[j * nx + 1];
		if (neumann && nbr[6] == -1) {
			f_comp[j * nx] = (-h_x * west - center + east) / (h_x * h_x);
		} else if (nbr[6] == -1) {
			f_comp[j * nx] = (2 * west - 3 * center + east) / (h_x * h_x);
		} else {
			f_comp[j * nx] = (west - 2 * center + east) / (h_x * h_x);
		}
	}
	// middle
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 0; j < ny; j++) {
			east   = u[j * nx + i - 1];
			center = u[j * nx + i];
			west   = u[j * nx + i + 1];

			f_comp[j * nx + i] = (west - 2 * center + east) / (h_x * h_x);
		}
	}
	// east
	for (int j = 0; j < ny; j++) {
		west   = u[j * nx + nx - 2];
		center = u[j * nx + nx - 1];
		east   = boundary_east[j];
		if (neumann && nbr[2] == -1) {
			f_comp[j * nx + nx - 1] = (west - center + h_x * east) / (h_x * h_x);
		} else if (nbr[2] == -1) {
			f_comp[j * nx + nx - 1] = (west - 3 * center + 2 * east) / (h_x * h_x);
		} else {
			f_comp[j * nx + nx - 1] = (west - 2 * center + east) / (h_x * h_x);
		}
	}
	// south
	for (int i = 0; i < nx; i++) {
		south  = boundary_south[i];
		center = u[i];
		north  = u[nx + i];
		if (neumann && nbr[4] == -1) {
			f_comp[i] += (-h_y * south - center + north) / (h_y * h_y);
		} else if (nbr[4] == -1) {
			f_comp[i] += (2 * south - 3 * center + north) / (h_y * h_y);
		} else {
			f_comp[i] += (south - 2 * center + north) / (h_y * h_y);
		}
	}
	// middle
	for (int i = 0; i < nx; i++) {
		for (int j = 1; j < ny - 1; j++) {
			south  = u[(j - 1) * nx + i];
			center = u[j * nx + i];
			north  = u[(j + 1) * nx + i];

			f_comp[j * nx + i] += (south - 2 * center + north) / (h_y * h_y);
		}
	}
	// north
	for (int i = 0; i < nx; i++) {
		south  = u[(ny - 2) * nx + i];
		center = u[(ny - 1) * nx + i];
		north  = boundary_north[i];
		if (neumann && nbr[0] == -1) {
			f_comp[(ny - 1) * nx + i] += (south - center + h_y * north) / (h_y * h_y);
		} else if (nbr[0] == -1) {
			f_comp[(ny - 1) * nx + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
		} else {
			f_comp[(ny - 1) * nx + i] += (south - 2 * center + north) / (h_y * h_y);
		}
	}
	resid = f - f_comp;
	return sqrt(pow(f - f_comp, 2).sum());
}

void Domain::solveWithInterface(const vector_type &gamma, vector_type &diff)
{
	auto gamma_view = gamma.getLocalView<Kokkos::HostSpace>();
	auto diff_view  = diff.getLocalView<Kokkos::HostSpace>();
	if (hasNbrNorth()) {
		boundary_north = valarray<double>(nx);
		int curr_i     = local_i[0];
		for (int i = 0; i < nx; i++) {
			boundary_north[i] = gamma_view(curr_i, 0);
			curr_i++;
		}
	}
	if (hasNbrEast()) {
		if (isRefinedEast()) {
			boundary_east = valarray<double>(ny);
			boundary_east_refined_right = valarray<double>(ny);
			int curr_i    = local_i[3];
			for (int i = 0; i < ny; i++) {
				boundary_east[i / 2] += gamma_view(curr_i, 0) / 2.0;
				boundary_east_refined_right[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
			curr_i = local_i[2];
			boundary_east_refined_left = valarray<double>(ny);
			for (int i = 0; i < ny; i++) {
				boundary_east[(i+nx) / 2 ] += gamma_view(curr_i, 0) / 2.0;
				boundary_east_refined_left[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
		} else {
			boundary_east = valarray<double>(ny);
			int curr_i    = local_i[2];
			for (int i = 0; i < ny; i++) {
				boundary_east[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
		}
	}
	if (hasNbrSouth()) {
		boundary_south = valarray<double>(nx);
		int curr_i     = local_i[4];
		for (int i = 0; i < nx; i++) {
			boundary_south[i] = gamma_view(curr_i, 0);
			curr_i++;
		}
	}
	if (hasNbrWest()) {
		boundary_west = valarray<double>(ny);
		int curr_i    = local_i[6];
		for (int i = 0; i < ny; i++) {
			boundary_west[i] = gamma_view(curr_i, 0);
			curr_i++;
		}
	}

	// solve
	solve();

	// if(has_east)cout <<"LOCAL Before\n";
	// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	// if(has_east)cout <<"LOCAL zero\n";
	// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	if (hasNbrNorth()) {
		int curr_i = local_i[0];
			valarray<double> diff   = getDiffNorth();
		for (int i = 0; i < nx; i++) {
			diff_view(curr_i, 0) += diff[i];
			curr_i++;
		}
	}

	if (hasNbrEast()) {
		if (isRefinedEast()) {
			valarray<double> diff   = getDiffEastRefinedLeft();
			int              curr_i = local_i[2];
			for (int i = 0; i < ny; i++) {
				diff_view(curr_i, 0) += diff[i];
				curr_i++;
			}
			diff   = getDiffEastRefinedRight();
			curr_i = local_i[3];
			for (int i = 0; i < ny; i++) {
				diff_view(curr_i, 0) += diff[i];
				curr_i++;
			}

            
		} else {
			int curr_i = local_i[2];
			valarray<double> diff   = getDiffEast();
			for (int i = 0; i < ny; i++) {
				diff_view(curr_i, 0) += diff[i];
				curr_i++;
			}
		}
	}
	if (hasNbrSouth()) {
		int              curr_i = local_i[4];
		valarray<double> diff   = getDiffSouth();
		for (int i = 0; i < nx; i++) {
			diff_view(curr_i, 0) += diff[i];
			curr_i++;
		}
	}
	if (hasNbrWest()) {
		int curr_i = local_i[6];
		valarray<double> diff   = getDiffWest();
		for (int i = 0; i < ny; i++) {
			diff_view(curr_i, 0) += diff[i];
			curr_i++;
		}
	}
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
valarray<double>                Domain::getDiffNorth()
{
	return boundary_north - valarray<double>(u[slice(nx * (nx - 1), nx, 1)]);
}
valarray<double> Domain::getDiffNorthRefinedLeft()
{
	return boundary_north - valarray<double>(u[slice(nx * (nx - 1), nx, 1)]);
}
valarray<double> Domain::getDiffNorthRefinedRight()
{
	return boundary_north - valarray<double>(u[slice(nx * (nx - 1), nx, 1)]);
}
valarray<double> Domain::getDiffEast()
{
	return boundary_east - valarray<double>(u[slice((nx - 1), nx, nx)]);
}
valarray<double> Domain::getDiffEastRefinedLeft()
{
	valarray<double> retval(nx);
	valarray<double> east = u[slice(nx-1,nx,nx)];
    int grid_i=nx-1;
    bool left = true;
	retval[nx - 1] = (east[nx - 1] + boundary_north[0]) / 2;
	for (int i = nx - 2; i >= 0; i--) {
		if (left) {
			retval[i] = (east[grid_i] + east[grid_i - 1]) / 2.0;
			left      = false;
		} else {
			retval[i] = (east[grid_i] + east[grid_i - 1]) / 2.0;
			left      = true;
			grid_i--;
		}
	}
	retval = 1.0 / 2.0 * (boundary_east_refined_left - retval);
	return retval;
}
valarray<double> Domain::getDiffEastRefinedRight()
{
	valarray<double> retval(nx);
	valarray<double> east = u[slice(nx - 1, nx, nx)];
	int              grid_i = 0;
	bool             left   = false;
	retval[0]               = (boundary_south[0] + east[0]) / 2;
	for (int i = 1; i < nx; i++) {
		if (left) {
			retval[i] = (east[grid_i] + east[grid_i + 1]) / 2.0;
			left      = false;
			grid_i++;
		} else {
			retval[i] = (east[grid_i] + east[grid_i + 1]) / 2.0;
			left      = true;
		}
	}
	retval = 1.0 / 2.0 * (boundary_east_refined_right - retval);
	return retval;
}
valarray<double> Domain::getDiffSouth()
{
	return boundary_south - valarray<double>(u[slice(0, nx, 1)]);
}
valarray<double> Domain::getDiffSouthRefinedLeft()
{

	valarray<double> retval(nx);
	/*
	retval[0]  = (boundary_west[0] + u[0]) / 2;
	int nleft  = (nx - 1) / 2 + (nx - 1) % 2;
	int nright = (nx - 1) / 2;
	retval[slice(1, nleft, 2)]
	= 3.0 / 4.0 * u[slice(0, nleft, 1)] + 1.0 / 4.0 * u[slice(1, nleft, 1)];
	retval[slice(2, nright, 2)]
	= 1.0 / 4.0 * u[slice(0, nright, 1)] + 3.0 / 4.0 * u[slice(1, nright, 1)];
	retval = 1.0 / 3.0 * retval;
    */
	return retval;
}
valarray<double> Domain::getDiffSouthRefinedRight()
{
	valarray<double> retval(nx);
    /*
	retval[nx-1] = (u[nx-1] + boundary_east[0]) / 2;
	int nleft  = (nx - 1) / 2;
	int nright = (nx - 1) / 2 + (nx - 1) % 2;
	retval[slice(nx-2*nright, nright, 2)]
	= 1.0 / 4.0 * u[slice(nx-2-nright, nright, 1)] + 3.0 / 4.0 * u[slice(nx-1-nright, nright, 1)];
	retval[slice(nx-1-2*nleft, nleft, 2)]
	= 3.0 / 4.0 * u[slice(nx-2-nleft, nleft, 1)] + 1.0 / 4.0 * u[slice(nx-1-nleft, nleft, 1)];
	retval = 1.0 / 3.0 * retval;
    */
	return retval;
}
valarray<double> Domain::getDiffWest()
{
	return boundary_west - valarray<double>(u[slice(0, nx, nx)]);
}
valarray<double> Domain::getDiffWestRefinedLeft()
{
	return boundary_west - valarray<double>(u[slice(0, nx, nx)]);
}
valarray<double> Domain::getDiffWestRefinedRight()
{
	return boundary_west - valarray<double>(u[slice(0, nx, nx)]);
}
