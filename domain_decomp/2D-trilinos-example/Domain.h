#ifndef DOMAIN_H
#define DOMAIN_H
#include "MyTypeDefs.h"
#include "DomainSignatureCollection.h"
#include <Teuchos_RCP.hpp>
#include <Tpetra_Import_decl.hpp>
#include <cmath>
#include <fftw3.h>
#include <valarray>
#include <array>
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
	std::array<int, 8> nbr           = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 8> local_i       = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 8> global_i      = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 8> iface_i       = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 8> iface_local_i = {-1, -1, -1, -1, -1, -1, -1, -1};
	Teuchos::RCP<map_type> domain_map;
	fftw_plan              plan1;
	fftw_plan              plan2;
	bool                   neumann = false;
	double                 x_start  = 0;
	double                 y_start  = 0;
	double                 x_length = 0;
	double                 y_length = 0;

	Domain() = default;
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
