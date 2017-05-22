#ifndef DOMAIN_H
#define DOMAIN_H
#include "DomainSignatureCollection.h"
#include "Side.h"
#include <array>
#include <cmath>
#include <valarray>
#include <HYPRE_sstruct_ls.h>
class Domain
{
	public:
	DomainSignature        ds;
	std::valarray<double>  f;
	std::valarray<double>  f_comp;
	std::valarray<double>  f_back;
	std::valarray<double>  resid;
	std::valarray<double>  u;
	std::valarray<double>  u_back;
	std::valarray<double>  exact;
	std::valarray<double>  tmp;
	std::valarray<double>  denom;
	std::valarray<double>  error;
	int                    n;
	double                 h_x;
	double                 h_y;
	std::valarray<double>  boundary_north;
	std::valarray<double>  boundary_south;
	std::valarray<double>  boundary_east;
	std::valarray<double>  boundary_west;
	std::valarray<double>  boundary_north_back;
	std::valarray<double>  boundary_south_back;
	std::valarray<double>  boundary_east_back;
	std::valarray<double>  boundary_west_back;
	std::array<int, 4> local_i        = {-1, -1, -1, -1};
	std::array<int, 4> local_i_center = {-1, -1, -1, -1};
	std::array<int, 4> local_i_left = {-1, -1, -1, -1};
	std::array<int, 4> local_i_right = {-1, -1, -1, -1};
	bool                   neumann = false;
	double                 x_start  = 0;
	double                 y_start  = 0;
	double                 x_length = 0;
	double                 y_length = 0;

	Domain() = default;
	Domain(DomainSignature ds, int n);
	~Domain();

	void setNeumann();
	double diffNorm();
	double diffNorm(double uavg, double eavg);
	double integrateU();
	double exactNorm();
	double fNorm();
	double exactNorm(double eavg);
	double integrateExact();
	double integrateF() { return f.sum() * h_x * h_y; }
	double residual();
	double sumBoundaryFlux();
	double integrateBoundaryFlux();
	double area();
	double integrateAU() { return f_comp.sum() * h_x * h_y; }
	void   swapResidSol();
	void   sumResidIntoSol();

	inline int &index(Side s)
	{
		return local_i[static_cast<int>(s)];
	}
	inline int &indexCenter(Side s) { return local_i_center[static_cast<int>(s)]; }
	inline int &indexRefinedLeft(Side s) { return local_i_left[static_cast<int>(s)]; }
	inline int &indexRefinedRight(Side s) { return local_i_right[static_cast<int>(s)]; }
	inline int &globalIndex(Side s) { return ds.index(s); }
	inline int &globalIndexCenter(Side s) { return ds.indexCenter(s); }
	inline int &globalIndexRefinedLeft(Side s) { return ds.indexRefinedLeft(s); }
	inline int &globalIndexRefinedRight(Side s) { return ds.indexRefinedRight(s); }
	inline int nbr(Side s) { return ds.nbr(s); }
	inline int nbrRight(Side      s) { return ds.nbrRight(s); }
	inline std::valarray<double> *getBoundaryPtr(Side s)
	{
		std::valarray<double> *retval = nullptr;
		switch (s) {
			case Side::north:
				retval = &boundary_north;
				break;
			case Side::east:
				retval = &boundary_east;
				break;
			case Side::south:
				retval = &boundary_south;
				break;
			case Side::west:
				retval = &boundary_west;
		}
		return retval;
	}
	inline std::valarray<double> getBoundary(const Side s) const
	{
		std::valarray<double> retval(n);
		switch (s) {
			case Side::north:
				retval = boundary_north;
				break;
			case Side::east:
				retval = boundary_east;
				break;
			case Side::south:
				retval = boundary_south;
				break;
			case Side::west:
				retval = boundary_west;
		}
		return retval;
	}
	inline bool hasNbr(Side s) const { return ds.hasNbr(s); }
	inline bool hasNonRefinedNbr(Side s) const
	{
		return ds.hasNbr(s) && !hasCoarseNbr(s) && !hasFineNbr(s);
	}
	inline bool hasFineNbr(Side s) const { return ds.hasFineNbr(s); }
	inline bool hasCoarseNbr(Side s) const { return ds.hasCoarseNbr(s); }
	inline bool isCoarseLeft(const Side s) const { return ds.leftOfCoarse(s); }
	inline bool leftToRight(Side s)
	{
		bool retval = (s == Side::north || s == Side::west);
		return retval;
	}
	void outputClaw(std::ostream &os);
	void setGridNbrs(HYPRE_SStructGrid &grid);
	void setStencils(HYPRE_SStructMatrix &A);
	void setAmrStencil(HYPRE_SStructGraph &graph);
	void setMatrixCoeffs(HYPRE_SStructMatrix &A);
	void fillRHS(HYPRE_SStructVector &b);
	void fillLHS(HYPRE_SStructVector &x);
	void fillExact(HYPRE_SStructVector &x);
	void saveLHS(HYPRE_SStructVector &x);
	void saveResid(HYPRE_SStructVector &r);
	void saveAU(HYPRE_SStructVector &b);
};
#endif
