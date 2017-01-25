#ifndef DOMAIN_H
#define DOMAIN_H
#include "MyTypeDefs.h"
#include <Epetra_Import.h>
#include <Teuchos_RCP.hpp>
#include <cmath>
#include <fftw3.h>
#include <valarray>
class Domain
{
	public:
	std::valarray<double>  f;
	std::valarray<double>  u;
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
	std::valarray<double>  orig_north;
	std::valarray<double>  orig_south;
	std::valarray<double>  orig_east;
	std::valarray<double>  orig_west;
	std::slice             north;
	std::slice             east;
	std::slice             south;
	std::slice             west;
	bool                   has_north = false;
	bool                   has_south = false;
	bool                   has_east  = false;
	bool                   has_west  = false;
	Teuchos::RCP<map_type> domain_map;
	fftw_plan              plan1;
	fftw_plan              plan2;

	Domain(std::valarray<double> &f, int nx, int ny, double h_x, double h_y)
	{
		this->f   = f;
		this->nx  = nx;
		this->ny  = ny;
		this->h_x = h_x;
		this->h_y = h_y;

		std::valarray<double> f_copy = f;
		tmp                          = std::valarray<double>(nx * ny);
		u                            = std::valarray<double>(nx * ny);
		denom                        = std::valarray<double>(nx * ny);
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
		plan1 = fftw_plan_r2r_2d(ny, nx, &f[0], &tmp[0], FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

		plan2 = fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);

		f = f_copy;
		// create slices
		north = std::slice(nx * (ny - 1), nx, 1);
		east  = std::slice((nx - 1), ny, nx);
		south = std::slice(0, nx, 1);
		west  = std::slice(0, ny, nx);

		orig_north = std::valarray<double>(nx);
		orig_north = f[north];
		orig_east  = std::valarray<double>(ny);
		orig_east  = f[east];
		orig_south = std::valarray<double>(nx);
		orig_south = f[south];
		orig_west  = std::valarray<double>(ny);
		orig_west  = f[west];
	}

	~Domain()
	{
		fftw_destroy_plan(plan1);
		fftw_destroy_plan(plan2);
	}

	void solve()
	{
		// f_copy = f;
		if (has_north) {
			f[north] = orig_north;
			f[north] -= 2 / (h_y * h_y) * boundary_north;
		}
		if (has_east) {
			f[east] = orig_east;
			f[east] -= 2 / (h_x * h_x) * boundary_east;
		}
		if (has_south) {
			f[south] = orig_south;
			f[south] -= 2 / (h_y * h_y) * boundary_south;
		}
		if (has_west) {
			f[west] = orig_west;
			f[west] -= 2 / (h_x * h_x) * boundary_west;
		}

		fftw_execute(plan1);

		tmp /= denom;

		fftw_execute(plan2);

		u /= 4 * nx * ny;
	}

	void solveWithInterface(const vector_type &gamma, vector_type &diff)
	{
		diff.PutScalar(0);
		Epetra_Import importer(*domain_map, gamma.Map());
		vector_type   local_vector(*domain_map, 1);
		local_vector.Import(gamma, importer, Insert);
		int curr_i = 0;
		if (has_north) {
			boundary_north = std::valarray<double>(&local_vector[0][curr_i], nx);
			curr_i += nx;
		}
		if (has_east) {
			boundary_east = std::valarray<double>(&local_vector[0][curr_i], ny);
			curr_i += ny;
		}
		if (has_south) {
			boundary_south = std::valarray<double>(&local_vector[0][curr_i], nx);
			curr_i += nx;
		}
		if (has_west) {
			boundary_west = std::valarray<double>(&local_vector[0][curr_i], ny);
		}

		// solve
		solve();

		local_vector.PutScalar(0.0);
		curr_i = 0;
		if (has_north) {
			for (int i = 0; i < nx; i++) {
				local_vector[0][curr_i] = u[nx * (ny - 1) + i];
				curr_i++;
			}
		}
		if (has_east) {
			for (int i = 0; i < ny; i++) {
				local_vector[0][curr_i] = u[(i + 1) * nx - 1];
				curr_i++;
			}
		}
		if (has_south) {
			for (int i = 0; i < nx; i++) {
				local_vector[0][curr_i] = u[i];
				curr_i++;
			}
		}
		if (has_west) {
			for (int i = 0; i < ny; i++) {
				local_vector[0][curr_i] = u[i * nx];
				curr_i++;
			}
		}
		diff.Export(local_vector, importer, Add);
		diff.Update(-2, gamma, 1);
		// diff.Print(std::cout);
		// sleep(1);
	}
};
#endif
