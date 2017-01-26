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
		if (has_north) {
			f_copy[std::slice(nx * (ny - 1), nx, 1)] -= 2 / (h_y * h_y) * boundary_north;
		}
		if (has_east) {
			f_copy[std::slice((nx - 1), ny, nx)] -= 2 / (h_x * h_x) * boundary_east;
		}
		if (has_south) {
			f_copy[std::slice(0, nx, 1)] -= 2 / (h_y * h_y) * boundary_south;
		}
		if (has_west) {
			f_copy[std::slice(0, ny, nx)] -= 2 / (h_x * h_x) * boundary_west;
		}

		fftw_execute(plan1);

		tmp /= denom;

		fftw_execute(plan2);

		u /= 4 * nx * ny;
	}

	void solveWithInterface(const vector_type &gamma, vector_type &diff)
	{
		// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
		// if(has_east)std::cout << "Gamma begin\n";
		// gamma.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		diff.update(1, gamma, 0);
		Tpetra::Import<> importer(diff.getMap(), domain_map);
		Tpetra::Export<> exporter(domain_map, diff.getMap());
		vector_type      local_vector(domain_map, 1);
		local_vector.doImport(diff, importer, Tpetra::CombineMode::INSERT);
		auto data_ptr = local_vector.getDataNonConst(0);
		int  curr_i   = 0;
		auto view     = local_vector.getLocalView<Kokkos::HostSpace>();
		if (has_north) {
			boundary_north = std::valarray<double>(nx);
			for (int i = 0; i < nx; i++) {
				boundary_north[i] = view(curr_i, 0);
				curr_i++;
			}
		}
		if (has_east) {
			boundary_east = std::valarray<double>(ny);
			for (int i = 0; i < ny; i++) {
				boundary_east[i] = view(curr_i, 0);
				curr_i++;
			}
		}
		if (has_south) {
			boundary_south = std::valarray<double>(nx);
			for (int i = 0; i < nx; i++) {
				boundary_south[i] = view(curr_i, 0);
				curr_i++;
			}
		}
		if (has_west) {
			boundary_west = std::valarray<double>(ny);
			for (int i = 0; i < ny; i++) {
				boundary_west[i] = view(curr_i, 0);
				curr_i++;
			}
		}

		// solve
		solve();

		// if(has_east)std::cout <<"LOCAL Before\n";
		// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		local_vector.scale(0);
		// if(has_east)std::cout <<"LOCAL zero\n";
		// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		curr_i = 0;
		if (has_north) {
			for (int i = 0; i < nx; i++) {
				view(curr_i, 0) = u[nx * (ny - 1) + i];
				curr_i++;
			}
		}
		if (has_east) {
			for (int i = 0; i < ny; i++) {
				view(curr_i, 0) = u[(i + 1) * nx - 1];
				curr_i++;
			}
		}
		if (has_south) {
			for (int i = 0; i < nx; i++) {
				view(curr_i, 0) = u[i];
				curr_i++;
			}
		}
		if (has_west) {
			for (int i = 0; i < ny; i++) {
				view(curr_i, 0) = u[i * nx];
				curr_i++;
			}
		}
		// if(has_east)std::cout <<"LOCAL AFEr\n";
		// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		diff.scale(0);
		// if(has_east)std::cout <<"DIFFBEFORE\n";
		// diff.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		diff.doExport(local_vector, importer, Tpetra::CombineMode::ADD);
		// if(has_east)std::cout <<"DIFF\n";
		// diff.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		// if(has_east)std::cout <<"GAMMA\n";
		// gamma.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
		diff.update(-2, gamma, 1);
	}
};
#endif
