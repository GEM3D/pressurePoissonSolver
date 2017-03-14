#ifndef DOMAIN_H
#define DOMAIN_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include <Teuchos_RCP.hpp>
#include <Tpetra_Import_decl.hpp>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <valarray>
enum class Dir { north, east, south, west };
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

	bool hasRefinedNbr(Dir dir);
	std::valarray<double> getDiffNorth();
	std::valarray<double> getDiffNorthRefinedLeft();
	std::valarray<double> getDiffNorthRefinedRight();
	std::valarray<double> getDiffEast();
	std::valarray<double> getDiffEastRefinedLeft();
	std::valarray<double> getDiffEastRefinedRight();
	std::valarray<double> getDiffSouth();
	std::valarray<double> getDiffSouthRefinedLeft();
	std::valarray<double> getDiffSouthRefinedRight();
	std::valarray<double> getDiffWest();
	std::valarray<double> getDiffWestRefinedLeft();
	std::valarray<double> getDiffWestRefinedRight();

	inline bool isRefinedNorth() { return nbr[1] != -1; }
	inline bool isRefinedEast() { return nbr[3] != -1; }
	inline bool isRefinedSouth() { return nbr[5] != -1; }
	inline bool isRefinedWest() { return nbr[7] != -1; }
	inline bool hasNbrNorth() { return nbr[0] != -1; }
	inline bool hasNbrEast() { return nbr[2] != -1; }
	inline bool hasNbrSouth() { return nbr[4] != -1; }
	inline bool hasNbrWest() { return nbr[6] != -1; }
};
#endif
