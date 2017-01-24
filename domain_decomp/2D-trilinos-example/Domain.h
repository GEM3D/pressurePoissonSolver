#ifndef DOMAIN_H
#define DOMAIN_H
#include "MyTypeDefs.h"
#include <Epetra_Import.h>
#include <Teuchos_RCP.hpp>
#include <valarray>
#include <cmath>
#include <fftw3.h>
class Domain
{
	public:
	std::valarray<double>  f;
	std::valarray<double>  f_copy;
	std::valarray<double>  u;
	std::valarray<double>  tmp;
	std::valarray<double>  denom;
	int                    nx;
	int                    ny;
	double                 h_x;
	double                 h_y;
	double *               boundary_north = nullptr;
	double *               boundary_south = nullptr;
	double *               boundary_east  = nullptr;
	double *               boundary_west  = nullptr;
	bool                   has_north      = false;
	bool                   has_south      = false;
	bool                   has_east       = false;
	bool                   has_west       = false;
	Teuchos::RCP<map_type> domain_map;

	Domain(std::valarray<double> &f, int nx, int ny, double h_x,
	       double h_y)
	{
		this->f   = f;
		this->nx  = nx;
		this->ny  = ny;
		this->h_x = h_x;
		this->h_y = h_y;

		f_copy = f;
        
		u     = f;
		tmp   = std::valarray<double>(nx * ny);
		denom = std::valarray<double>(nx * ny);
		for (int yi = 1; yi <= ny; yi++) {
			for (int xi = 1; xi <= nx; xi++) {
				denom[nx * (yi - 1) + xi - 1]
				= -4 / (h_x * h_x) * std::pow(std::sin(xi * M_PI / (2 * nx)), 2)
				  - 4 / (h_y * h_y) * std::pow(std::sin(yi * M_PI / (2 * ny)), 2);
			}
		}
	}

	void solve()
	{
		for (int y_i = 0; y_i < ny; y_i++) {
			for (int x_i = 0; x_i < nx; x_i++) {
				f_copy[y_i * nx + x_i] = f[y_i * nx + x_i];
				if (y_i == ny - 1 && has_north) {
					f_copy[y_i * nx + x_i] += -2.0 / (h_y * h_y) * boundary_north[x_i];
				}
				if (x_i == nx - 1 && has_east) {
					f_copy[y_i * nx + x_i] += -2.0 / (h_x * h_x) * boundary_east[y_i];
				}
				if (y_i == 0 && has_south) {
					f_copy[y_i * nx + x_i] += -2.0 / (h_y * h_y) * boundary_south[x_i];
				}
				if (x_i == 0 && has_west) {
					f_copy[y_i * nx + x_i] += -2.0 / (h_x * h_x) * boundary_west[y_i];
				}
			}
		}

		// create fftw plans
		fftw_plan plan1 = fftw_plan_r2r_2d(ny, nx, &f_copy[0], &tmp[0], FFTW_RODFT10,
		                                   FFTW_RODFT10, FFTW_MEASURE);

		fftw_plan plan2
		= fftw_plan_r2r_2d(ny, nx, &tmp[0], &u[0], FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);

		fftw_execute(plan1);

		tmp /= denom;

		fftw_execute(plan2);

		u /= 4 * nx * ny;

		fftw_destroy_plan(plan1);
		fftw_destroy_plan(plan2);
	}

	void solveWithInterface(const vector_type &gamma, vector_type &diff)
	{
		diff.PutScalar(0);
		Epetra_Import importer(*domain_map, gamma.Map());
		vector_type   local_vector(*domain_map, 1);
		local_vector.Import(gamma, importer, Insert);
		int curr_i = 0;
		if (has_north) {
			boundary_north = &local_vector[0][curr_i];
			curr_i += nx;
		}
		if (has_east) {
			boundary_east = &local_vector[0][curr_i];
			curr_i += ny;
		}
		if (has_south) {
			boundary_south = &local_vector[0][curr_i];
			curr_i += nx;
		}
		if (has_west) {
			boundary_west = &local_vector[0][curr_i];
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
