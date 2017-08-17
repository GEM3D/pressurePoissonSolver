#ifndef DOMAIN_H
#define DOMAIN_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include "Side.h"
#include "Solver.h"
#include <Teuchos_RCP.hpp>
#include <Tpetra_Import_decl.hpp>
#include <array>
#include <bitset>
#include <cmath>
#include <valarray>
class Domain
{
	public:
	Solver *              solver;
	DomainSignature       ds;
	std::valarray<double> f;
	std::valarray<double> f_comp;
	std::valarray<double> f_back;
	std::valarray<double> resid;
	std::valarray<double> u;
	std::valarray<double> u_back;
	std::valarray<double> exact;
	std::valarray<double> tmp;
	std::valarray<double> denom;
	std::valarray<double> error;
	int                   n;
	double                h_x;
	double                h_y;
	std::valarray<double> boundary_north;
	std::valarray<double> boundary_south;
	std::valarray<double> boundary_east;
	std::valarray<double> boundary_west;
	std::valarray<double> boundary_north_back;
	std::valarray<double> boundary_south_back;
	std::valarray<double> boundary_east_back;
	std::valarray<double> boundary_west_back;
	Teuchos::RCP<map_type> domain_map;
	bool                   neumann    = false;
	double                 x_start    = 0;
	double                 y_start    = 0;
	double                 x_length   = 0;
	double                 y_length   = 0;

	Domain() = default;
	Domain(DomainSignature ds, int n);
	~Domain();

	void plan();
	void planDirichlet();
	void planNeumann();
	void solve();
	void putGhostCells(vector_type &ghost);
	double residual(vector_type &ghost);
	double residual();
	void solveWithInterface(const vector_type &gamma);
	void getDiff(vector_type &diff);
	double diffNorm();
	double diffNorm(double uavg, double eavg);
	double integrateU();
	double exactNorm();
	double fNorm();
	double exactNorm(double eavg);
	double integrateExact();
	double integrateF() { return f.sum() * h_x * h_y; }
	double sumBoundaryFlux();
	double integrateBoundaryFlux();
	double area();
	double integrateAU() { return f_comp.sum() * h_x * h_y; }
	void   swapResidSol();
	void   sumResidIntoSol();

	std::valarray<double> getDiff(const Side s) const;
	std::valarray<double> getDiffCombined(const Side s) const;
	std::valarray<double> getDiffCombinedLeft(const Side s) const;
	std::valarray<double> getDiffCombinedRight(const Side s) const;
	std::valarray<double> getDiffFine(const Side s) const;
	std::valarray<double> getDiffFineToCoarse(const Side s) const;
	std::valarray<double> getDiffFineToCoarseLeft(const Side s) const;
	std::valarray<double> getDiffFineToCoarseRight(const Side s) const;
	std::valarray<double> getDiffCoarse(const Side s) const;
	std::valarray<double> getDiffCoarseToFineLeft(const Side s) const;
	std::valarray<double> getDiffCoarseToFineRight(const Side s) const;

	std::valarray<double> getSide(const Side s) const;
	std::valarray<double> getInnerSide(const Side s) const;
	std::valarray<double> getSideFine(const Side s) const;
	std::valarray<double> getSideFineLeft(const Side s) const;
	std::valarray<double> getSideFineRight(const Side s) const;
	std::valarray<double> getInnerSideFine(const Side s) const;
	std::valarray<double> getInnerSideFineLeft(const Side s) const;
	std::valarray<double> getInnerSideFineRight(const Side s) const;
	std::valarray<double> getSideCoarse(const Side s) const;
	std::valarray<double> getSideCoarseCombined(const Side s) const;
	inline std::valarray<double> getSideCoarseLeft(const Side s) const
	{
		if (s == Side::west || s == Side::north) {
			return getSideCoarseRelativeLeft(s);
		} else {
			return getSideCoarseRelativeRight(s);
		}
	}
	inline std::valarray<double> getSideCoarseRight(const Side s) const
	{
		if (s == Side::west || s == Side::north) {
			return getSideCoarseRelativeRight(s);
		} else {
			return getSideCoarseRelativeLeft(s);
		}
	}

	std::valarray<double> getSideCoarseRelativeLeft(const Side s) const;
	std::valarray<double> getSideCoarseRelativeRight(const Side s) const;
	void fillBoundary(Side s, const single_vector_type &gamma);
	void fillDiffVector(Side s, single_vector_type &diff);
	void fillGhostVector(Side s, single_vector_type &diff);
	void getFluxDiff(vector_type &flux);
	void putCoords(vector_type &xy);
	void fillCoords(Side s,vector_type &xy);
	void fillFluxVector(Side s, single_vector_type &diff);

	inline int &index(Side s) { return ds.index(s); }
	inline int &indexCenter(Side s) { return ds.indexCenter(s); }
	inline int &indexRefinedLeft(Side s) { return ds.indexRefinedLeft(s); }
	inline int &indexRefinedRight(Side s) { return ds.indexRefinedRight(s); }
	inline int &globalIndex(Side s) { return ds.globalIndex(s); }
	inline int &globalIndexCenter(Side s) { return ds.globalIndexCenter(s); }
	inline int &globalIndexRefinedLeft(Side s) { return ds.globalIndexRefinedLeft(s); }
	inline int &globalIndexRefinedRight(Side s) { return ds.globalIndexRefinedRight(s); }
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
	inline bool hasFineNbr(Side s) const { return ds.hasFineNbr(s); }
	inline bool hasCoarseNbr(Side s) const { return ds.hasCoarseNbr(s); }
	inline bool isCoarseLeft(const Side s) const { return ds.leftOfCoarse(s); }
	inline bool leftToRight(Side s)
	{
		bool retval = (s == Side::north || s == Side::west);
		return retval;
	}
	void outputClaw(std::ostream &os);
};
#endif
