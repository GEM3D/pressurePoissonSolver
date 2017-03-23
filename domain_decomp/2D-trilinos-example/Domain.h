#ifndef DOMAIN_H
#define DOMAIN_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include "Side.h"
#include <Teuchos_RCP.hpp>
#include <Tpetra_Import_decl.hpp>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <valarray>
class Domain
{
	public:
	DomainSignature        ds;
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
	std::valarray<double>  boundary_east_refined_left;
	std::valarray<double>  boundary_east_refined_right;
	std::valarray<double>  boundary_west;
	std::array<int, 8>  nbr           = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 12> local_i       = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 12> global_i      = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 12> iface_i       = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 12> iface_local_i = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
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

	std::valarray<double> getStencil(Side s, Tilt t = Tilt::center);
	std::valarray<double> getSide(Side s);
	std::valarray<double> getSideCoarseLeft(Side s);
	std::valarray<double> getSideCoarseRight(Side s);
	std::valarray<double> getSideFine(Side s);
	void fillBoundary(Side s, const single_vector_type &gamma);
	void fillDiffVector(Side s, single_vector_type &diff, bool weight = true);

	inline bool hasCoarseNbr(Side s)
	{
		bool retval;
		switch (s) {
			case Side::north:
				retval = ds.nbr_coarse[0];
				break;
			case Side::east:
				retval = ds.nbr_coarse[1];
				break;
			case Side::south:
				retval = ds.nbr_coarse[2];
				break;
			case Side::west:
				retval = ds.nbr_coarse[3];
		}
		return retval;
	}
	inline bool hasNbr(Side s)
	{
		bool retval;
		switch (s) {
			case Side::north:
				retval = nbr[0] != -1;
				break;
			case Side::east:
				retval = nbr[2] != -1;
				break;
			case Side::south:
				retval = nbr[4] != -1;
				break;
			case Side::west:
				retval = nbr[6] != -1;
		}
		return retval;
	}
	inline bool hasFineNbr(Side s)
	{
		bool retval;
		switch (s) {
			case Side::north:
				retval = nbr[1] != -1;
				break;
			case Side::east:
				retval = nbr[3] != -1;
				break;
			case Side::south:
				retval = nbr[5] != -1;
				break;
			case Side::west:
				retval = nbr[7] != -1;
		}
		return retval;
	}
};
#endif
