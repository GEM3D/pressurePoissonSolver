#ifndef DOMAIN_H
#define DOMAIN_H
#include "MyTypeDefs.h"
#include "DomainSignatureCollection.h"
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
	std::valarray<double>  resid;
	std::valarray<double>  u;
	std::valarray<double>  exact;
	std::valarray<double>  tmp;
	std::valarray<double>  denom;
	std::valarray<double>  error;
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
	double                 x_start  = 0;
	double                 y_start  = 0;
	double                 x_length = 0;
	double                 y_length = 0;

	Domain() {}
	Domain(DomainSignature ds, int nx, int ny, double h_x, double h_y);
	~Domain();

	void planDirichlet();
	void planNeumann();
	void solve();
	void putGhostCells(vector_type &ghost);
	double residual(vector_type &ghost);
	void solveWithInterface(const vector_type &gamma, vector_type &diff);
	double diffNorm();
	double diffNorm(double uavg, double eavg);
	double uSum();
	double exactNorm();
	double fNorm();
	double exactNorm(double eavg);
	double exactSum();
};
#endif
