#include "Domain.h"
#include "StencilHelper.h"
#include "CoeffHelper.h"
#include <array>
#include <ostream>
#include <valarray>
#include <vector>
//#include "FishpackSolver.h"
using namespace std;
Domain::Domain(DomainSignature ds, int n)
{
	this->n  = n;
	this->h_x = ds.x_length / n;
	this->h_y = ds.y_length / n;
	this->ds  = ds;

#if NDEBUG
	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "h:           " << this->h_x << endl;
	cerr << "I start at:  " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "Length:     " << ds.x_length << ", " << ds.y_length << "\n";
	cerr << "North: " << ds.nbr(Side::north) << ", " << ds.nbrRight(Side::north) << "\n";
	cerr << "Idx:   " << ds.index(Side::north) << ", " << ds.indexCenter(Side::north) << "\n";
	cerr << "       " << indexRefinedLeft(Side::north) << ", "
	     << indexRefinedRight(Side::north) << "\n";
	cerr << "East:  " << ds.nbr(Side::east) << ", " << ds.nbrRight(Side::east) << "\n";
	cerr << "Idx:   " << ds.index(Side::east) << ", " << ds.indexCenter(Side::east) << "\n";
	cerr << "       " << indexRefinedLeft(Side::east) << ", "
	     << indexRefinedRight(Side::east) << "\n";
	cerr << "South: " << ds.nbr(Side::south) << ", " << ds.nbrRight(Side::south) << "\n";
	cerr << "Idx:   " << ds.index(Side::south) << ", " << ds.indexCenter(Side::south) << "\n";
	cerr << "       " << indexRefinedLeft(Side::south) << ", "
	     << indexRefinedRight(Side::south) << "\n";
	cerr << "West:  " << ds.nbr(Side::west) << ", " << ds.nbrRight(Side::west) << "\n";
	cerr << "Idx:   " << ds.index(Side::west) << ", " << ds.indexCenter(Side::west) << "\n";
	cerr << "       " << indexRefinedLeft(Side::west) << ", "
	     << indexRefinedRight(Side::west) << "\n";
	cerr << "\n";
#endif
	f      = valarray<double>(n * n);
	f_back = valarray<double>(n * n);
	f_comp = valarray<double>(n * n);
	exact  = valarray<double>(n * n);
	u      = valarray<double>(n * n);
	u_back = valarray<double>(n * n);
	denom  = valarray<double>(n * n);
	resid  = valarray<double>(n * n);

	boundary_north = valarray<double>(n);
	boundary_south = valarray<double>(n);
	boundary_east  = valarray<double>(n);
	boundary_west  = valarray<double>(n);

	boundary_north_back = valarray<double>(n);
	boundary_south_back = valarray<double>(n);
	boundary_east_back  = valarray<double>(n);
	boundary_west_back  = valarray<double>(n);

	x_start  = ds.x_start;
	y_start  = ds.y_start;
	x_length = ds.x_length;
	y_length = ds.y_length;
}

Domain::~Domain()
{
#if NDEBUG
	cerr << "I am Domain: " << ds.id << "\n";
	cerr << "h:           " << this->h_x << endl;
	cerr << "I start at:  " << ds.x_start << ", " << ds.y_start << "\n";
	cerr << "Length:     " << ds.x_length << ", " << ds.y_length << "\n";
	cerr << "North: " << ds.nbr(Side::north) << ", " << ds.nbrRight(Side::north) << "\n";
	cerr << "Idx:   " << index(Side::north) << ", " << indexCenter(Side::north) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::north) << ", " << hasCoarseNbr(Side::north) << ", "
	     << isCoarseLeft(Side::north) << endl;
	cerr << "East:  " << ds.nbr(Side::east) << ", " << ds.nbrRight(Side::east) << "\n";
	cerr << "Idx:   " << index(Side::east) << ", " << indexCenter(Side::east) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::east) << ", " << hasCoarseNbr(Side::east) << ", "
	     << isCoarseLeft(Side::east) << endl;
	cerr << "South: " << ds.nbr(Side::south) << ", " << ds.nbrRight(Side::south) << "\n";
	cerr << "Idx:   " << index(Side::south) << ", " << indexCenter(Side::south) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::south) << ", " << hasCoarseNbr(Side::south) << ", "
	     << isCoarseLeft(Side::south) << endl;
	cerr << "West:  " << ds.nbr(Side::west) << ", " << ds.nbrRight(Side::west) << "\n";
	cerr << "Idx:   " << index(Side::west) << ", " << indexCenter(Side::west) << "\n";
	cerr << "Flags: " << hasFineNbr(Side::west) << ", " << hasCoarseNbr(Side::west) << ", "
	     << isCoarseLeft(Side::west) << endl;
	cerr << "\n";
#endif
}

void Domain::setNeumann()
{
	neumann = true;
}

double Domain::diffNorm()
{
	error = exact - u;
	return sqrt(pow(exact - u, 2).sum());
}
double Domain::diffNorm(double uavg, double eavg)
{
	error = exact - u - eavg + uavg;
	return sqrt(pow(exact - u - eavg + uavg, 2).sum());
}
double Domain::integrateU() { return u.sum() * h_x * h_y; }
double Domain::exactNorm() { return sqrt((exact * exact).sum()); }
double Domain::fNorm() { return sqrt((f * f).sum()); }
double Domain::residual() { return sqrt((resid*resid).sum()); }
double Domain::exactNorm(double eavg) { return sqrt(pow(exact - eavg, 2).sum()); }
double                          Domain::integrateExact() { return exact.sum()*h_x*h_y; }
double Domain::area() { return x_length * y_length; }
double Domain::integrateBoundaryFlux()
{
	double sum = 0;
	if (!hasNbr(Side::north)) {
		sum += boundary_north.sum() * h_x;
	}
	if (!hasNbr(Side::east)) {
		sum -= boundary_east.sum() * h_y;
	}
	if (!hasNbr(Side::south)) {
		sum -= boundary_south.sum() * h_x;
	}
	if (!hasNbr(Side::west)) {
		sum += boundary_west.sum() * h_y;
	}
    return sum;
}
void Domain::outputClaw(std::ostream &os)
{
	const string tab = "\t";
	os << ds.id << tab << "grid_number" << endl;
	os << ds.refine_level << tab << "AMR_level" << endl;
	os << 0 << tab << "block_number" << endl;
	os << 0 << tab << "mpi_rank" << endl;
	os << n << tab << "mx" << endl;
	os << n << tab << "my" << endl;
	os << x_start << tab << "xlow" << endl;
	os << y_start << tab << "ylow" << endl;
	os << h_x << tab << "dx" << endl;
	os << h_y << tab << "dy" << endl;
	os << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int loc = j + i * n;
			os << u[loc] << tab << error[loc] << tab << resid[loc] * h_x * h_y << endl;
		}
		os << endl;
	}
}
void Domain::setGridNbrs(HYPRE_SStructGrid &grid) {
	int index_map[2] = {0, 1};
	int index_dir[2] = {1, 1};
	if (hasNonRefinedNbr(Side::north)) {
		int ilower[2]     = {0, n};
		int iupper[2]     = {n - 1, n};
		int nbr_ilower[2] = {0, 0};
		int nbr_iupper[2] = {n - 1, 0};
		HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::north), nbr_ilower,
		                                 nbr_iupper, index_map, index_dir);
	}
	if (hasNonRefinedNbr(Side::east)) {
		int ilower[2]     = {n, 0};
		int iupper[2]     = {n, n - 1};
		int nbr_ilower[2] = {0, 0};
		int nbr_iupper[2] = {0, n - 1};
		HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::east), nbr_ilower,
		                                 nbr_iupper, index_map, index_dir);
	}
	if (hasNonRefinedNbr(Side::south)) {
		int ilower[2]     = {0, -1};
		int iupper[2]     = {n - 1, -1};
		int nbr_ilower[2] = {0, n - 1};
		int nbr_iupper[2] = {n - 1, n - 1};
		HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::south), nbr_ilower,
		                                 nbr_iupper, index_map, index_dir);
	}
	if (hasNonRefinedNbr(Side::west)) {
		int ilower[2]     = {-1, 0};
		int iupper[2]     = {-1, n - 1};
		int nbr_ilower[2] = {n - 1, 0};
		int nbr_iupper[2] = {n - 1, n - 1};
		HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::west), nbr_ilower,
		                                 nbr_iupper, index_map, index_dir);
	}
}

void Domain::setAmrStencil(HYPRE_SStructGraph &graph) {
	auto getCoarseHelper = [&](Side s) {
		CoarseStencilHelper *ret = nullptr;
        switch(s){
			case Side::north:
				ret = new NorthCSH(this);
				break;
			case Side::east:
				ret = new EastCSH(this);
				break;
			case Side::south:
				ret = new SouthCSH(this);
				break;
			case Side::west:
				ret = new WestCSH(this);
        }
		return ret;
	};
	auto getFineHelper = [&](Side s) {
		FineStencilHelper *ret = nullptr;
		switch (s) {
			case Side::north:
				ret = new NorthFSH(this);
				break;
			case Side::east:
				ret = new EastFSH(this);
				break;
			case Side::south:
				ret = new SouthFSH(this);
				break;
			case Side::west:
				ret = new WestFSH(this);
		}
		return ret;
	};

	Side s = Side::north;
	do {
		if (hasCoarseNbr(s)) {
			CoarseStencilHelper *sh = getCoarseHelper(s);
			while (!sh->done()) {
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexCoarseLeft, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexCoarseCenter, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexCoarseRight, 0);
				++(*sh);
			}
			delete sh;
		}

		if (hasFineNbr(s)) {
			FineStencilHelper *sh = getFineHelper(s);
			HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, ds.id,
			                             sh->nbrIndexFineOtherLeft, 0);
			while (!sh->done()) {
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexFineLeftIn, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexFineLeftOut, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexFineRightIn, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, sh->nbr(),
				                             sh->nbrIndexFineRightOut, 0);
				++(*sh);
			}
			HYPRE_SStructGraphAddEntries(graph, ds.id, sh->currIndex, 0, ds.id,
			                             sh->nbrIndexFineOtherRight, 0);
			delete sh;
		}
		s++;
	} while (s != Side::north);
}
void Domain::setMatrixCoeffs(HYPRE_SStructMatrix &A)
{
	// middle
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		int stencil_indeces[3] = {0, 3, 4};
		vector<array<double, 3>> coeffs;
		array<double, 3>         this_coeffs;
		this_coeffs[0] = -2.0 / (h_y * h_y);
		this_coeffs[1] = 1.0 / (h_y * h_y);
		this_coeffs[2] = 1.0 / (h_y * h_y);
		for (int yi = 1; yi < n - 1; yi++) {
			for (int xi = 0; xi < n; xi++) {
				coeffs.push_back(this_coeffs);
			}
		}
		int lower[2] = {0, 1};
		int upper[2] = {n - 1, n - 2};
		HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
		                                  &coeffs[0][0]);
	}
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		int stencil_indeces[3] = {0, 1, 2};
		vector<array<double, 3>> coeffs;
		array<double, 3>         this_coeffs;
		this_coeffs[0] = -2.0 / (h_x * h_x);
		this_coeffs[1] = 1.0 / (h_x * h_x);
		this_coeffs[2] = 1.0 / (h_x * h_x);
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 1; xi < n - 1; xi++) {
				coeffs.push_back(this_coeffs);
			}
		}
		int lower[2] = {1, 0};
		int upper[2] = {n - 2, n - 1};
		HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
		                                  &coeffs[0][0]);
	}

    //patch boundaries
	Side s = Side::north;
	do {
		CoeffHelper ch = CoeffHelper(this, s);
		while (!ch.done()) {
			HYPRE_SStructMatrixAddToValues(A, ds.id, ch.curr_index, 0, ch.size, ch.indices,
			                               ch.values);
			HYPRE_SStructMatrixAddToValues(A, ds.id, ch.curr_index, 0, ch.off_size, ch.off_indices,
			                               ch.off_values);
			++ch;
		}
		s++;
	} while (s != Side::north);
	if (ds.id == 0 && neumann) {
		int    idx[2]     = {1, 1};
		int    indices[5] = {0, 1, 2, 3, 4};
		double values[5]  = {1, 0, 0, 0, 0};
		HYPRE_SStructMatrixSetValues(A, ds.id, idx, 0, 5, indices, values);
	}
}
void Domain::fillRHS(HYPRE_SStructVector &b)
{
	valarray<double> f_copy = f;
	if (!hasNbr(Side::north)) {
		if (neumann) {
			f_copy[slice(n * (n - 1), n, 1)] -= 1 / h_y * boundary_north;
		} else {
			f_copy[slice(n * (n - 1), n, 1)] -= 2 / (h_y * h_y) * boundary_north;
		}
	}
	if (!hasNbr(Side::east)) {
		if (neumann) {
			f_copy[slice((n - 1), n, n)] -= 1 / h_x * boundary_east;
		} else {
			f_copy[slice((n - 1), n, n)] -= 2 / (h_x * h_x) * boundary_east;
		}
	}
	if (!hasNbr(Side::south)) {
		if (neumann) {
			f_copy[slice(0, n, 1)] += 1 / h_y * boundary_south;
		} else {
			f_copy[slice(0, n, 1)] -= 2 / (h_y * h_y) * boundary_south;
		}
	}
	if (!hasNbr(Side::west)) {
		if (neumann) {
			f_copy[slice(0, n, n)] += 1 / h_x * boundary_west;
		} else {
			f_copy[slice(0, n, n)] -= 2 / (h_x * h_x) * boundary_west;
		}
	}
    int lower[2]={0,0};
    int upper[2]={n-1,n-1};
	HYPRE_SStructVectorSetBoxValues(b, ds.id, lower, upper, 0, &f_copy[0]);
}
void Domain::fillLHS(HYPRE_SStructVector &x) {
	int lower[2] = {0, 0};
	int upper[2] = {n-1, n-1};
	HYPRE_SStructVectorSetBoxValues(x, ds.id, lower, upper, 0, &u[0]);
}
void Domain::fillExact(HYPRE_SStructVector &x) {
	int lower[2] = {0, 0};
	int upper[2] = {n-1, n-1};
	HYPRE_SStructVectorSetBoxValues(x, ds.id, lower, upper, 0, &exact[0]);
}
void Domain::saveLHS(HYPRE_SStructVector &x)
{
	int lower[2] = {0, 0};
	int upper[2] = {n-1, n-1};
	HYPRE_SStructVectorGetBoxValues(x, ds.id, lower, upper, 0, &u[0]);
}
void Domain::saveResid(HYPRE_SStructVector &r)
{
	int lower[2] = {0, 0};
	int upper[2] = {n-1, n-1};
	HYPRE_SStructVectorGetBoxValues(r, ds.id, lower, upper, 0, &resid[0]);
}
void Domain::saveAU(HYPRE_SStructVector &b)
{
	int lower[2] = {0, 0};
	int upper[2] = {n-1, n-1};
	HYPRE_SStructVectorGetBoxValues(b, ds.id, lower, upper, 0, &f_comp[0]);
	if (!hasNbr(Side::north)) {
		if (neumann) {
			f_comp[slice(n * (n - 1), n, 1)] += 1 / h_y * boundary_north;
		} else {
			f_comp[slice(n * (n - 1), n, 1)] += 2 / (h_y * h_y) * boundary_north;
		}
	}
	if (!hasNbr(Side::east)) {
		if (neumann) {
			f_comp[slice((n - 1), n, n)] += 1 / h_x * boundary_east;
		} else {
			f_comp[slice((n - 1), n, n)] += 2 / (h_x * h_x) * boundary_east;
		}
	}
	if (!hasNbr(Side::south)) {
		if (neumann) {
			f_comp[slice(0, n, 1)] -= 1 / h_y * boundary_south;
		} else {
			f_comp[slice(0, n, 1)] += 2 / (h_y * h_y) * boundary_south;
		}
	}
	if (!hasNbr(Side::west)) {
		if (neumann) {
			f_comp[slice(0, n, n)] -= 1 / h_x * boundary_west;
		} else {
			f_comp[slice(0, n, n)] += 2 / (h_x * h_x) * boundary_west;
		}
	}
}
