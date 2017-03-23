#include "Domain.h"
using namespace std;
Domain::Domain(DomainSignature ds, int nx, int ny, double h_x, double h_y)
{
	this->nx  = nx;
	this->ny  = ny;
	this->h_x = h_x / ds.refine_level;
	this->h_y = h_y / ds.refine_level;

	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "h:           " << this->h_x << endl;
	cerr << "I start at:  " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "Length:     " << ds.x_length << ", " << ds.y_length << "\n";
	cerr << "North: " << ds.nbr[0] << ", " << ds.nbr[1] << "\n";
	cerr << "Idx:   " << ds.global_i[0] << ", " << ds.global_i[1] << ", " << ds.global_i[2] << "\n";
	cerr << "East:  " << ds.nbr[2] << ", " << ds.nbr[3] << "\n";
	cerr << "Idx:   " << ds.global_i[3] << ", " << ds.global_i[4] << ", " << ds.global_i[5] << "\n";
	cerr << "South: " << ds.nbr[4] << ", " << ds.nbr[5] << "\n";
	cerr << "Idx:   " << ds.global_i[6] << ", " << ds.global_i[7] << ", " << ds.global_i[8] << "\n";
	cerr << "West:  " << ds.nbr[6] << ", " << ds.nbr[7] << "\n";
	cerr << "Idx:   " << ds.global_i[9] << ", " << ds.global_i[10] << ", " << ds.global_i[11]
	     << "\n";
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

	for (int q = 0; q < 12; q++) {
		if (ds.global_i[q] != -1) {
			global_i[q] = ds.global_i[q] * nx;
			iface_i[q]  = ds.global_i[q] * 22;
		}
	}

	nbr = ds.nbr;

	x_start  = ds.x_start;
	y_start  = ds.y_start;
	x_length = ds.x_length;
	y_length = ds.y_length;

	this->ds = ds;
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
		denom[slice(yi * ny, nx, 1)] = -4 / (h_x * h_x) * pow(sin((yi + 1) * M_PI / (2 * nx)), 2);
	}

	valarray<double> ones(ny);
	ones = 1;
	for (int xi = 0; xi < ny; xi++) {
		denom[slice(xi, ny, nx)]
		-= 4 / (h_y * h_y) * pow(sin((xi + 1) * M_PI / (2 * ny)), 2) * ones;
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
	plan1 = fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], y_transform, x_transform, FFTW_MEASURE);

	plan2
	= fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], y_transform_inv, x_transform_inv, FFTW_MEASURE);

	// create denom vector
	if (!hasNbr(Side::north) && !hasNbr(Side::south)) {
		for (int yi = 0; yi < nx; yi++) {
			denom[slice(yi * ny, nx, 1)] = -4 / (h_x * h_x) * pow(sin(yi * M_PI / (2 * nx)), 2);
		}
	} else if (!hasNbr(Side::south) || !hasNbr(Side::north)) {
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

	if (!hasNbr(Side::east) && !hasNbr(Side::west)) {
		for (int xi = 0; xi < ny; xi++) {
			denom[slice(xi, ny, nx)] -= 4 / (h_y * h_y) * pow(sin(xi * M_PI / (2 * ny)), 2) * ones;
		}
	} else if (!hasNbr(Side::west) || !hasNbr(Side::east)) {
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
	if (!hasNbr(Side::north) && neumann) {
		f_copy[slice(nx * (ny - 1), nx, 1)] -= 1 / h_y * boundary_north;
	} else {
		f_copy[slice(nx * (ny - 1), nx, 1)] -= 2 / (h_y * h_y) * boundary_north;
	}
	if (!hasNbr(Side::east) && neumann) {
		f_copy[slice((nx - 1), ny, nx)] -= 1 / h_x * boundary_east;
	} else {
		f_copy[slice((nx - 1), ny, nx)] -= 2 / (h_x * h_x) * boundary_east;
	}
	if (!hasNbr(Side::south) && neumann) {
		f_copy[slice(0, nx, 1)] += 1 / h_y * boundary_south;
	} else {
		f_copy[slice(0, nx, 1)] -= 2 / (h_y * h_y) * boundary_south;
	}
	if (!hasNbr(Side::west) && neumann) {
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
	auto left_ptr  = ghost.getVectorNonConst(0);
	auto right_ptr = ghost.getVectorNonConst(1);
	if (hasNbr(Side::north)) {
		fillDiffVector(Side::north, *right_ptr, true);
	}
	if (hasNbr(Side::east)) {
		fillDiffVector(Side::east, *right_ptr, true);
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

    /*
	if (hasNbr(Side::north)) {
		fillBoundary(Side::north, *left_ptr);
	}
	if (hasNbr(Side::east)) {
		fillBoundary(Side::east, *left_ptr);
	}
	if (hasNbr(Side::south)) {
		fillBoundary(Side::south, *right_ptr);
	}
	if (hasNbr(Side::west)) {
		fillBoundary(Side::west, *right_ptr);
	}
    */

	valarray<double> f_comp = valarray<double>(nx * ny);
	// integrate in x secton
	double center, north, east, south, west;
	// west
	for (int j = 0; j < ny; j++) {
		west   = boundary_west[j];
		center = u[j * nx];
		east   = u[j * nx + 1];
		if (neumann && !hasNbr(Side::west)) {
			f_comp[j * nx] = (-h_x * west - center + east) / (h_x * h_x);
		} else {
			f_comp[j * nx] = (2 * west - 3 * center + east) / (h_x * h_x);
        }
        /*
		} else if (!hasNbr(Side::west)) {
		} else if (hasFineNbr(Side::west)) {
			// TODO
			// f_comp[j * nx] = (1.0 / 3.0 * west - 4.0 / 3.0 * center + east) / (h_x * h_x);
		} else if (hasCoarseNbr(Side::west)) {
			// TODO
			f_comp[j * nx] = (2.6118 * west - 3.6105 * center + east) / (h_x * h_x);
            cerr << west << ", " << center << ", " << east << ", " << f[j*nx]*h_x*h_x << ", " <<f_comp[j*nx]*h_x*h_x << endl;
			f_comp[j * nx] = (west - 2 * center + east) / (h_x * h_x);
		}
        */
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
		if (neumann && !hasNbr(Side::east)) {
			f_comp[j * nx + nx - 1] = (west - center + h_x * east) / (h_x * h_x);
		} else {
			f_comp[j * nx + nx - 1] = (west - 3 * center + 2 * east) / (h_x * h_x);
        }
        /*
		} else if (!hasNbr(Side::east)) {
		} else if (hasFineNbr(Side::east)) {
			f_comp[j * nx + nx - 1] = (west - 7.0 / 3.0 * center + 4.0 / 3.0 * east) / (h_x * h_x);
		} else if (hasCoarseNbr(Side::east)) {
			// TODO
			// f_comp[j * nx + nx - 1] = (west - 4.0 / 3.0 * center + 1.0 / 3.0 * east) / (h_x *
			// h_x);
			f_comp[j * nx + nx - 1] = (west - 2 * center + east) / (h_x * h_x);
		}
        */
	}
	// south
	for (int i = 0; i < nx; i++) {
		south  = boundary_south[i];
		center = u[i];
		north  = u[nx + i];
		if (neumann && !hasNbr(Side::south)) {
			f_comp[i] += (-h_y * south - center + north) / (h_y * h_y);
		} else {
			f_comp[i] += (2 * south - 3 * center + north) / (h_y * h_y);
        }
        /*
		} else if (!hasNbr(Side::south)) {
			f_comp[i] += (south - 2 * center + north) / (h_y * h_y);
		}
        */
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
		if (neumann && !hasNbr(Side::north)) {
			f_comp[(ny - 1) * nx + i] += (south - center + h_y * north) / (h_y * h_y);
		} else {
			f_comp[(ny - 1) * nx + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
        }
        /*
		} else if (!hasNbr(Side::north)) {
			f_comp[(ny - 1) * nx + i] += (south - 2 * center + north) / (h_y * h_y);
		}
        */
	}
	if (hasFineNbr(Side::east)) {
		for (int j = 0; j < ny; j++) {
            f_comp[j*nx+nx-1]=f[j*nx+nx-1];
		}
	}
	resid = f - f_comp;
	return sqrt(pow(f - f_comp, 2).sum());
}

void Domain::solveWithInterface(const vector_type &gamma, vector_type &diff)
{
	auto vec_ptr  = gamma.getVector(0);
	auto diff_ptr = diff.getVectorNonConst(0);

	if (hasNbr(Side::north)) {
		fillBoundary(Side::north, *vec_ptr);
	}
	if (hasNbr(Side::east)) {
		fillBoundary(Side::east, *vec_ptr);
	}
	if (hasNbr(Side::south)) {
		fillBoundary(Side::south, *vec_ptr);
	}
	if (hasNbr(Side::west)) {
		fillBoundary(Side::west, *vec_ptr);
	}

	// solve
	solve();

	if (hasNbr(Side::north)) {
		fillDiffVector(Side::north, *diff_ptr);
	}
	if (hasNbr(Side::east)) {
		fillDiffVector(Side::east, *diff_ptr);
	}
	if (hasNbr(Side::south)) {
		fillDiffVector(Side::south, *diff_ptr);
	}
	if (hasNbr(Side::west)) {
		fillDiffVector(Side::west, *diff_ptr);
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
valarray<double> Domain::getSide(Side s){
	valarray<double> retval(nx);
	switch (s) {
		case Side::north:
			retval    = u[slice(nx * (nx - 1), nx, 1)];
			break;
		case Side::east:
			retval    = u[slice(nx - 1, nx, nx)];
			break;
		case Side::south:
			retval    = u[slice(0, nx, 1)];
			break;
		case Side::west:
			retval    = u[slice(0, nx, nx)];
	}
    return retval;
}
valarray<double> Domain::getSideCoarseLeft(Side s)
{
	valarray<double> retval(nx);
	valarray<double> side = getSide(s);
	switch (s) {
		case Side::north:
			retval[0] = (boundary_west[nx - 1] + side[0])/2.0;
			for (int i = 1; i < nx; i++) {
                if(i%2==1){
					retval[i] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
			break;
		case Side::east:
			retval[nx-1] = (boundary_north[nx-1] + side[nx-1])/2.0;
			for (int i = nx; i < 2*nx-1; i++) {
                if(i%2==1){
					retval[i-nx] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i-nx] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
			break;
		case Side::south:
			retval[nx-1] = (boundary_east[0] + side[nx-1])/2.0;
			for (int i = nx; i < 2*nx-1; i++) {
                if(i%2==1){
					retval[i-nx] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i-nx] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
			break;
		case Side::west:
			retval[0] = (boundary_south[0] + side[0])/2.0;
			for (int i = 1; i < nx; i++) {
                if(i%2==1){
					retval[i] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
	}
    return retval;
}
valarray<double> Domain::getSideCoarseRight(Side s)
{
	valarray<double> retval(nx);
	valarray<double> side = getSide(s);
	switch (s) {
		case Side::north:
			retval[nx-1] = (boundary_east[nx - 1] + side[nx-1])/2.0;
			for (int i = nx; i < 2*nx-1; i++) {
                if(i%2==1){
					retval[i-nx] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i-nx] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
			break;
		case Side::east:
			retval[0] = (boundary_south[nx-1] + side[0])/2.0;
			for (int i = 1; i < nx; i++) {
                if(i%2==1){
					retval[i] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
			break;
		case Side::south:
			retval[0] = (boundary_west[0] + side[0])/2.0;
			for (int i = 1; i < nx; i++) {
                if(i%2==1){
					retval[i] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
			break;
		case Side::west:
			retval[nx-1] = (boundary_north[0] + side[nx-1])/2.0;
			for (int i = nx; i < 2*nx-1; i++) {
                if(i%2==1){
					retval[i-nx] = (3.0 * side[(i - 1) / 2] + side[(i - 1) / 2 + 1]) / 4.0;
				}else{
					retval[i-nx] = (side[(i - 1) / 2] + 3.0 * side[(i - 1) / 2 + 1]) / 4.0;
				}
			}
	}
    return retval;
}
valarray<double> Domain::getSideFine(Side s)
{
	valarray<double> retval(nx);
	valarray<double> side = getSide(s);
	switch (s) {
		case Side::north:
			// left case
			if (global_i[1] != -1) {
				for (int i = 0; i < nx; i++) {
					retval[i / 2] += side[i];
				}
			}else{
				for (int i = 0; i < nx; i++) {
					retval[nx - 1 - i / 2] += side[nx - 1 - i];
				}
            }
			break;
		case Side::east:
			// left case
			if (global_i[4] != -1) {
				for (int i = 0; i < nx; i++) {
					retval[nx - 1 - i / 2] += side[nx - 1 - i];
				}
			} else {
				for (int i = 0; i < nx; i++) {
					retval[i / 2] += side[i];
				}
			}
			break;
		case Side::south:
			// left case
			if (global_i[7] != -1) {
				for (int i = 0; i < nx; i++) {
					retval[nx - 1 - i / 2] += side[nx - 1 - i];
				}
			} else {
				for (int i = 0; i < nx; i++) {
					retval[i / 2] += side[i];
				}
			}
			break;
		case Side::west:
			// left case
			if (global_i[10] != -1) {
				for (int i = 0; i < nx; i++) {
					retval[nx - 1 - i / 2] += side[nx - 1 - i];
				}
			} else {
				for (int i = 0; i < nx; i++) {
					retval[i / 2] += side[i];
				}
			}
	}
    retval/=2;
    return retval;
}
valarray<double> Domain::getStencil(Side s,Tilt t)
{
		return getSide(s);
}

void Domain::fillBoundary(Side s, const single_vector_type &gamma)
{
	int               curr_i;
	valarray<double> *boundary_ptr;
    if(hasCoarseNbr(s)){
		switch (s) {
			case Side::north:
				if (global_i[1] != -1) {
					curr_i = global_i[1];
				} else {
					curr_i = global_i[2];
				}
				boundary_ptr = &boundary_north;
				break;
			case Side::east:
				if (global_i[4] != -1) {
					curr_i = global_i[4];
				} else {
					curr_i = global_i[5];
				}
				boundary_ptr = &boundary_east;
				break;
			case Side::south:
				if (global_i[7] != -1) {
					curr_i = global_i[7];
				} else {
					curr_i = global_i[8];
				}
				boundary_ptr = &boundary_south;
				break;
			case Side::west:
				if (global_i[10] != -1) {
					curr_i = global_i[10];
				} else {
					curr_i = global_i[11];
				}
				boundary_ptr = &boundary_west;
		}
	} else {
		switch (s) {
			case Side::north:
				curr_i       = global_i[0];
				boundary_ptr = &boundary_north;
				break;
			case Side::east:
				curr_i       = global_i[3];
				boundary_ptr = &boundary_east;
				break;
			case Side::south:
				curr_i       = global_i[6];
				boundary_ptr = &boundary_south;
				break;
			case Side::west:
				curr_i       = global_i[9];
				boundary_ptr = &boundary_west;
		}
	}
	valarray<double> &boundary   = *boundary_ptr;
	auto              gamma_view = gamma.getLocalView<Kokkos::HostSpace>();

	for (int i = 0; i < nx; i++) {
		boundary[i] = gamma_view(curr_i, 0);
		curr_i++;
	}
}
void Domain::fillDiffVector(Side s, single_vector_type &diff, bool residual)
{
	// weight=false;
	auto diff_view = diff.getLocalView<Kokkos::HostSpace>();
	if (hasFineNbr(s)) {
		valarray<double> center = getSide(s);
		valarray<double> left   = getSideCoarseLeft(s);
		valarray<double> right  = getSideCoarseRight(s);
		if (!residual) {
			center *= 2.0 / 3.0;
			left *= 2.0 / 3.0;
			right *= 2.0 / 3.0;
		}
		int curr_i;
		int curr_i_left;
		int curr_i_right;
		switch (s) {
			case Side::north:
				curr_i       = global_i[0];
				curr_i_left  = global_i[1];
				curr_i_right = global_i[2];
				break;
			case Side::east:
				curr_i       = global_i[3];
				curr_i_left  = global_i[4];
				curr_i_right = global_i[5];
				break;
			case Side::south:
				curr_i       = global_i[6];
				curr_i_left  = global_i[7];
				curr_i_right = global_i[8];
				break;
			case Side::west:
				curr_i       = global_i[9];
				curr_i_left  = global_i[10];
				curr_i_right = global_i[11];
		}
		for (int i = 0; i < nx; i++) {
			diff_view(curr_i, 0) += center[i];
			if (residual) {
				diff_view(curr_i_left, 0) += left[i];
				diff_view(curr_i_right, 0) += right[i];
			}
			curr_i++;
			curr_i_left++;
			curr_i_right++;
		}
    }else if(hasCoarseNbr(s)){
		valarray<double> center = getSideFine(s);
		valarray<double> side = getSide(s);
        if(!residual){
			center *= 4.0 / 3.0;
			side *= 2.0;
		}
		int curr_i_center;
		int curr_i;
		switch (s) {
			case Side::north:
				curr_i_center = global_i[0];
				if (global_i[1] != -1) {
					curr_i = global_i[1];
				} else {
					curr_i = global_i[2];
				}
				break;
			case Side::east:
                curr_i_center = global_i[3];
				if (global_i[4] != -1) {
					curr_i = global_i[4];
				} else {
					curr_i = global_i[5];
				}
				break;
			case Side::south:
                curr_i_center = global_i[6];
				if (global_i[7] != -1) {
					curr_i = global_i[7];
				} else {
					curr_i = global_i[8];
				}
                break;
			case Side::west:
                curr_i_center = global_i[9];
				if (global_i[10] != -1) {
					curr_i = global_i[10];
				} else {
					curr_i = global_i[11];
				}
		}
		for (int i = 0; i < nx; i++) {
			diff_view(curr_i_center, 0) += center[i];
			diff_view(curr_i, 0) += side[i];
			curr_i++;
			curr_i_center++;
		}
	} else {
		valarray<double> side = getSide(s);
		int curr_i;
		switch (s) {
			case Side::north:
				curr_i = global_i[0];
				break;
			case Side::east:
				curr_i = global_i[3];
				break;
			case Side::south:
				curr_i = global_i[6];
				break;
			case Side::west:
				curr_i = global_i[9];
		}
		for (int i = 0; i < nx; i++) {
			diff_view(curr_i, 0) += side[i];
			curr_i++;
		}
	}
}
