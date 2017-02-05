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
	int                    nbr_north = -1;
	int                    nbr_east  = -1;
	int                    nbr_south = -1;
	int                    nbr_west  = -1;
	int                    local_i_north;
	int                    local_i_east;
	int                    local_i_south;
	int                    local_i_west;
	int                    global_i_north;
	int                    global_i_east;
	int                    global_i_south;
	int                    global_i_west;
	int                    iface_i_north;
	int                    iface_i_east;
	int                    iface_i_south;
	int                    iface_i_west;
	int                    iface_local_i_north;
	int                    iface_local_i_east;
	int                    iface_local_i_south;
	int                    iface_local_i_west;
	Teuchos::RCP<map_type> domain_map;
	fftw_plan              plan1;
	fftw_plan              plan2;

	Domain() {}
	Domain(std::valarray<double> f, std::valarray<double> exact, int nx, int ny, double h_x,
	       double h_y)
	{
		this->f     = f;
		this->exact = exact;
		this->nx    = nx;
		this->ny    = ny;
		this->h_x   = h_x;
		this->h_y   = h_y;

		f_copy = f;
		tmp    = std::valarray<double>(nx * ny);
		u      = std::valarray<double>(nx * ny);
		denom  = std::valarray<double>(nx * ny);
		for (int xi = 1; xi <= nx; xi++) {
			denom[std::slice((xi - 1) * ny, nx, 1)]
			= -4 / (h_x * h_x) * std::pow(std::sin(xi * M_PI / (2 * nx)), 2);
		}
		std::valarray<double> ones(ny);
		ones = 1;
		for (int yi = 1; yi <= ny; yi++) {
			denom[std::slice((yi - 1), ny, nx)]
			-= 4 / (h_y * h_y) * std::pow(std::sin(yi * M_PI / (2 * ny)), 2) * ones;
		}
		// create fftw plans
		plan1
		= fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

		plan2 = fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);
	}

	~Domain()
	{
		fftw_destroy_plan(plan1);
		fftw_destroy_plan(plan2);
	}

	void solve()
	{
		f_copy = f;
		if (nbr_north != -1) {
			f_copy[std::slice(nx * (ny - 1), nx, 1)] -= 2 / (h_y * h_y) * boundary_north;
		}
		if (nbr_east != -1) {
			f_copy[std::slice((nx - 1), ny, nx)] -= 2 / (h_x * h_x) * boundary_east;
		}
		if (nbr_south != -1) {
			f_copy[std::slice(0, nx, 1)] -= 2 / (h_y * h_y) * boundary_south;
		}
		if (nbr_west != -1) {
			f_copy[std::slice(0, ny, nx)] -= 2 / (h_x * h_x) * boundary_west;
		}

		fftw_execute(plan1);

		tmp /= denom;

		fftw_execute(plan2);

		u /= 4 * nx * ny;
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
	double exactNorm() { return std::sqrt((exact * exact).sum()); }
};
#endif
