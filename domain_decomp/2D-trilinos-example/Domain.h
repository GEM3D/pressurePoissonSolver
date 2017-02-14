#ifndef DOMAIN_H
#define DOMAIN_H
#include "MyTypeDefs.h"
#include <Teuchos_RCP.hpp>
#include <Tpetra_Import_decl.hpp>
#include <cmath>
#include <fftw3.h>
#include <valarray>
class Domain
{
	public:
	std::valarray<double>  f;
	std::valarray<double>  f_copy;
	std::valarray<double>  u;
	std::valarray<double>  exact;
	std::valarray<double>  tmp;
	std::valarray<double>  denom;
	int                    nx;
	int                    ny;
	double                 h_x;
	double                 h_y;
	std::valarray<double>  boundary_north;
	std::valarray<double>  boundary_south;
	std::valarray<double>  boundary_east;
	std::valarray<double>  boundary_west;
	int                    nbr_north           = -1;
	int                    nbr_east            = -1;
	int                    nbr_south           = -1;
	int                    nbr_west            = -1;
	int                    local_i_north       = -1;
	int                    local_i_east        = -1;
	int                    local_i_south       = -1;
	int                    local_i_west        = -1;
	int                    global_i_north      = -1;
	int                    global_i_east       = -1;
	int                    global_i_south      = -1;
	int                    global_i_west       = -1;
	int                    iface_i_north       = -1;
	int                    iface_i_east        = -1;
	int                    iface_i_south       = -1;
	int                    iface_i_west        = -1;
	int                    iface_local_i_north = -1;
	int                    iface_local_i_east  = -1;
	int                    iface_local_i_south = -1;
	int                    iface_local_i_west  = -1;
	Teuchos::RCP<map_type> domain_map;
	fftw_plan              plan1;
	fftw_plan              plan2;
	bool                   neumann = false;

	Domain() {}
	Domain(int nx, int ny, double h_x, double h_y)
	{
		this->nx  = nx;
		this->ny  = ny;
		this->h_x = h_x;
		this->h_y = h_y;

		f      = std::valarray<double>(nx * ny);
		f_copy = std::valarray<double>(nx * ny);
		exact  = std::valarray<double>(nx * ny);
		tmp    = std::valarray<double>(nx * ny);
		u      = std::valarray<double>(nx * ny);
		denom  = std::valarray<double>(nx * ny);

		boundary_north = std::valarray<double>(nx);
		boundary_south = std::valarray<double>(nx);
		boundary_east  = std::valarray<double>(ny);
		boundary_west  = std::valarray<double>(ny);
	}

	~Domain()
	{
		fftw_destroy_plan(plan1);
		fftw_destroy_plan(plan2);
	}

	void planDirichlet()
	{
        //create plan
		plan1
		= fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

		plan2 = fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);
        //create denom vector

		for (int yi = 0; yi < nx; yi++) {
			denom[std::slice(yi * ny, nx, 1)]
			= -4 / (h_x * h_x) * std::pow(std::sin((yi + 1) * M_PI / (2 * nx)), 2);
		}

		std::valarray<double> ones(ny);
		ones = 1;
		for (int xi = 0; xi < ny; xi++) {
			denom[std::slice(xi, ny, nx)]
			-= 4 / (h_y * h_y) * std::pow(std::sin((xi + 1) * M_PI / (2 * ny)), 2) * ones;
		}
	}

	void planNeumann()
	{
		neumann                       = true;
		fftw_r2r_kind x_transform     = FFTW_RODFT10;
		fftw_r2r_kind x_transform_inv = FFTW_RODFT01;
		fftw_r2r_kind y_transform     = FFTW_RODFT10;
		fftw_r2r_kind y_transform_inv = FFTW_RODFT01;
		if (nbr_east == -1 && nbr_west == -1) {
			x_transform     = FFTW_REDFT10;
			x_transform_inv = FFTW_REDFT01;
		} else if (nbr_west == -1) {
			x_transform     = FFTW_REDFT11;
			x_transform_inv = FFTW_REDFT11;
		} else if (nbr_east == -1) {
			x_transform     = FFTW_RODFT11;
			x_transform_inv = FFTW_RODFT11;
		}
		if (nbr_north == -1 && nbr_south == -1) {
			y_transform     = FFTW_REDFT10;
			y_transform_inv = FFTW_REDFT01;
		} else if (nbr_south == -1) {
			y_transform     = FFTW_REDFT11;
			y_transform_inv = FFTW_REDFT11;
		} else if (nbr_north == -1) {
			y_transform     = FFTW_RODFT11;
			y_transform_inv = FFTW_RODFT11;
		}
		plan1
		= fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], y_transform, x_transform, FFTW_MEASURE);

		plan2
		= fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], y_transform_inv, x_transform_inv, FFTW_MEASURE);

		// create denom vector
		if (nbr_north == -1 && nbr_south == -1) {
			for (int yi = 0; yi < nx; yi++) {
				denom[std::slice(yi * ny, nx, 1)]
				= -4 / (h_x * h_x) * std::pow(std::sin(yi * M_PI / (2 * nx)), 2);
			}
		} else if (nbr_south == -1 || nbr_north == -1) {
			for (int yi = 0; yi < nx; yi++) {
				denom[std::slice(yi * ny, nx, 1)]
				= -4 / (h_x * h_x) * std::pow(std::sin((yi + 0.5) * M_PI / (2 * nx)), 2);
			}
		} else {
			for (int yi = 0; yi < nx; yi++) {
				denom[std::slice(yi * ny, nx, 1)]
				= -4 / (h_x * h_x) * std::pow(std::sin((yi + 1) * M_PI / (2 * nx)), 2);
			}
		}

		std::valarray<double> ones(ny);
		ones = 1;

		if (nbr_east == -1 && nbr_west == -1) {
			for (int xi = 0; xi < ny; xi++) {
				denom[std::slice(xi, ny, nx)]
				-= 4 / (h_y * h_y) * std::pow(std::sin(xi * M_PI / (2 * ny)), 2) * ones;
			}
		} else if (nbr_west == -1 || nbr_east == -1) {
			for (int xi = 0; xi < ny; xi++) {
				denom[std::slice(xi, ny, nx)]
				-= 4 / (h_y * h_y) * std::pow(std::sin((xi + 0.5) * M_PI / (2 * ny)), 2) * ones;
			}
		} else {
			for (int xi = 0; xi < ny; xi++) {
				denom[std::slice(xi, ny, nx)]
				-= 4 / (h_y * h_y) * std::pow(std::sin((xi + 1) * M_PI / (2 * ny)), 2) * ones;
			}
		}
	}

	void solve()
	{
		f_copy = f;
        if(nbr_north==-1&&neumann){
			f_copy[std::slice(nx * (ny - 1), nx, 1)] -= 2 / (h_y * h_y) * boundary_north;
        }else{
			f_copy[std::slice(nx * (ny - 1), nx, 1)] -= 2 / (h_y * h_y) * boundary_north;
        }
		if (nbr_east == -1&&neumann) {
			f_copy[std::slice((nx - 1), ny, nx)] -= 2 / (h_x * h_x) * boundary_east;
		}else{
			f_copy[std::slice((nx - 1), ny, nx)] -= 2 / (h_x * h_x) * boundary_east;
        }
		if (nbr_south == -1&&neumann) {
			f_copy[std::slice(0, nx, 1)] -= 2 / (h_y * h_y) * boundary_south;
		}else{
			f_copy[std::slice(0, nx, 1)] -= 2 / (h_y * h_y) * boundary_south;
        }
		if (nbr_west == -1&&neumann) {
			f_copy[std::slice(0, ny, nx)] -= 2 / (h_x * h_x) * boundary_west;
		}else{
			f_copy[std::slice(0, ny, nx)] -= 2 / (h_x * h_x) * boundary_west;
        }

		fftw_execute(plan1);

		tmp /= denom;

		if (neumann && nbr_north == -1 && nbr_east == -1 && nbr_south == -1 && nbr_west == -1) {
			tmp[0] = 0;
		}

		fftw_execute(plan2);

		u /= 4 * nx * ny;
	}

	void putGhostCells(vector_type &ghost)
	{
		auto ghost_view = ghost.getLocalView<Kokkos::HostSpace>();
		if (nbr_north != -1) {
			int curr_i = local_i_north;
			for (int i = 0; i < nx; i++) {
				ghost_view(curr_i, 1) += u[nx * (ny - 1) + i];
				curr_i++;
			}
		}
		if (nbr_east != -1) {
			int curr_i = local_i_east;
			for (int i = 0; i < ny; i++) {
				ghost_view(curr_i, 1) += u[(i + 1) * nx - 1];
				curr_i++;
			}
		}
		if (nbr_south != -1) {
			int curr_i = local_i_south;
			for (int i = 0; i < nx; i++) {
				ghost_view(curr_i, 0) += u[i];
				curr_i++;
			}
		}
		if (nbr_west != -1) {
			int curr_i = local_i_west;
			for (int i = 0; i < ny; i++) {
				ghost_view(curr_i, 0) += u[i * nx];
				curr_i++;
			}
		}
	}

	double residual(vector_type &ghost)
	{
		auto ghost_view = ghost.getLocalView<Kokkos::HostSpace>();
		if (nbr_north != -1) {
			boundary_north = std::valarray<double>(nx);
			int curr_i     = local_i_north;
			for (int i = 0; i < nx; i++) {
				boundary_north[i] = ghost_view(curr_i, 0);
				curr_i++;
			}
		}
		if (nbr_east != -1) {
			boundary_east = std::valarray<double>(ny);
			int curr_i    = local_i_east;
			for (int i = 0; i < ny; i++) {
				boundary_east[i] = ghost_view(curr_i, 0);
				curr_i++;
			}
		}
		if (nbr_south != -1) {
			boundary_south = std::valarray<double>(nx);
			int curr_i     = local_i_south;
			for (int i = 0; i < nx; i++) {
				boundary_south[i] = ghost_view(curr_i, 1);
				curr_i++;
			}
		}
		if (nbr_west != -1) {
			boundary_west = std::valarray<double>(ny);
			int curr_i    = local_i_west;
			for (int i = 0; i < ny; i++) {
				boundary_west[i] = ghost_view(curr_i, 1);
				curr_i++;
			}
		}

		std::valarray<double> f_comp = std::valarray<double>(nx * ny);
		// integrate in x directon
		double center, north, east, south, west;
		// west
		for (int j = 0; j < ny; j++) {
			west   = boundary_west[j];
			center = u[j * nx];
			east   = u[j * nx + 1];
			if (neumann && nbr_west == -1) {
				f_comp[j * nx] = (-h_x * west - center + east) / (h_x * h_x);
			} else if (nbr_west == -1) {
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
			if (neumann && nbr_east == -1) {
				f_comp[j * nx + nx - 1] = (west - center + h_x * east) / (h_x * h_x);
			} else if (nbr_east == -1) {
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
			if (neumann && nbr_south == -1) {
				f_comp[i] += (-h_y * south - center + north) / (h_y * h_y);
			} else if (nbr_south == -1) {
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
			if (neumann && nbr_north == -1) {
				f_comp[(ny - 1) * nx + i] += (south - center + h_y * north) / (h_y * h_y);
			} else if (nbr_north == -1) {
				f_comp[(ny - 1) * nx + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
			} else {
				f_comp[(ny - 1) * nx + i] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
		return sqrt(pow(f - f_comp, 2).sum());
	}

	void solveWithInterface(const vector_type &gamma, vector_type &diff)
	{
		auto gamma_view = gamma.getLocalView<Kokkos::HostSpace>();
		auto diff_view  = diff.getLocalView<Kokkos::HostSpace>();
		if (nbr_north != -1) {
			boundary_north = std::valarray<double>(nx);
			int curr_i     = local_i_north;
			for (int i = 0; i < nx; i++) {
				boundary_north[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
		}
		if (nbr_east != -1) {
			boundary_east = std::valarray<double>(ny);
			int curr_i    = local_i_east;
			for (int i = 0; i < ny; i++) {
				boundary_east[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
		}
		if (nbr_south != -1) {
			boundary_south = std::valarray<double>(nx);
			int curr_i     = local_i_south;
			for (int i = 0; i < nx; i++) {
				boundary_south[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
		}
		if (nbr_west != -1) {
			boundary_west = std::valarray<double>(ny);
			int curr_i    = local_i_west;
			for (int i = 0; i < ny; i++) {
				boundary_west[i] = gamma_view(curr_i, 0);
				curr_i++;
			}
		}

		// solve
		solve();

		// if(has_east)std::cout <<"LOCAL Before\n";
		// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		// if(has_east)std::cout <<"LOCAL zero\n";
		// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		if (nbr_north != -1) {
			int curr_i = local_i_north;
			for (int i = 0; i < nx; i++) {
				diff_view(curr_i, 0) += u[nx * (ny - 1) + i];
				curr_i++;
			}
		}
		if (nbr_east != -1) {
			int curr_i = local_i_east;
			for (int i = 0; i < ny; i++) {
				diff_view(curr_i, 0) += u[(i + 1) * nx - 1];
				curr_i++;
			}
		}
		if (nbr_south != -1) {
			int curr_i = local_i_south;
			for (int i = 0; i < nx; i++) {
				diff_view(curr_i, 0) += u[i];
				curr_i++;
			}
		}
		if (nbr_west != -1) {
			int curr_i = local_i_west;
			for (int i = 0; i < ny; i++) {
				diff_view(curr_i, 0) += u[i * nx];
				curr_i++;
			}
		}
	}
	double diffNorm() { return std::sqrt(std::pow(exact - u, 2).sum()); }
	double diffNorm(double uavg, double eavg)
	{
		return std::sqrt(std::pow(exact - u - eavg + uavg, 2).sum());
	}
	double uSum() { return u.sum(); }
	double exactNorm() { return std::sqrt((exact * exact).sum()); }
	double fNorm() { return std::sqrt((f * f).sum()); }
	double exactNorm(double eavg) { return std::sqrt(std::pow(exact - eavg, 2).sum()); }
	double                  exactSum() { return exact.sum(); }
};
#endif
