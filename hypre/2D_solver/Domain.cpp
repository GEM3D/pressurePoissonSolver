#include "Domain.h"
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
	if (hasCoarseNbr(Side::west)) {
        Side s = Side::west;
        if(isCoarseLeft(s)){
			// first stencil is special
			int curr_idx[2] = {0, n - 1};
			int nbr_idx1[2] = {n - 1, n-3};
			int nbr_idx2[2] = {n - 1, n-2};
			int nbr_idx3[2] = {n - 1, n-1};
			// right
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx3, 0);
			curr_idx[1] = n - 2;
			// left
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx3, 0);
			// set rest of entries
			for (int i = 0; i < n - 2; i++) {
				curr_idx[1] = i;
				nbr_idx1[1] = n / 2 + i / 2 - 1;
				nbr_idx2[1] = n / 2 + i / 2;
				nbr_idx3[1] = n / 2 + i / 2 + 1;
				HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx3, 0);

			}
        }else{
			// first stencil is special
			int curr_idx[2] = {0, 0};
			int nbr_idx1[2] = {n - 1, 0};
			int nbr_idx2[2] = {n - 1, 1};
			int nbr_idx3[2] = {n - 1, 2};
			// left
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx3, 0);
			curr_idx[1] = 1;
			// right
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx3, 0);
			// set rest of entries
			for (int i = 2; i < n; i++) {
				curr_idx[1] = i;
				nbr_idx1[1] = i / 2 - 1;
				nbr_idx2[1] = i / 2;
				nbr_idx3[1] = i / 2 + 1;
				HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
				HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx3, 0);

			}
		}
	}
	if (hasFineNbr(Side::east)) {
		Side      s           = Side::east;
		int       curr_idx[2] = {n - 1, n - 1};
		int       nbr_idx1[2] = {0, n - 1};
		int       nbr_idx2[2] = {1, n - 1};
		const int start_i     = n - 1;
		for (int i = 0; i < n; i++) {
			curr_idx[1] = start_i - i / 2;
			nbr_idx1[1] = start_i - i;
			nbr_idx2[1] = start_i - i;
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx1, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbr(s), nbr_idx2, 0);
		}
		for (int i = 0; i < n; i++) {
			curr_idx[1] = start_i - n / 2 - i / 2;
			nbr_idx1[1] = start_i - i;
			nbr_idx2[1] = start_i - i;
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbrRight(s), nbr_idx1, 0);
			HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, nbrRight(s), nbr_idx2, 0);
		}
		int this_idx[2] = {n - 1, n - 3};
		curr_idx[1]     = n - 1;
		HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, ds.id, this_idx, 0);
		this_idx[1] = 2;
		curr_idx[1] = 0;
		HYPRE_SStructGraphAddEntries(graph, ds.id, curr_idx, 0, ds.id, this_idx, 0);
	}
}
void Domain::setMatrixCoeffs(HYPRE_SStructMatrix &A) {
	// south
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		int stencil_indeces[3] = {0, 3, 4};
		vector<array<double, 3>> coeffs;
		array<double, 3>         this_coeffs;
		if (hasNbr(Side::south)) {
			this_coeffs[0] = -2.0 / (h_y * h_y);
			this_coeffs[1] = 1.0 / (h_y * h_y);
			this_coeffs[2] = 1.0 / (h_y * h_y);
		} else {
			this_coeffs[0] = -3.0 / (h_y * h_y);
			this_coeffs[1] = 0;
			this_coeffs[2] = 1.0 / (h_y * h_y);
		}
		for (int xi = 0; xi < n; xi++) {
			coeffs.push_back(this_coeffs);
		}
		int lower[2] = {0, 0};
		int upper[2] = {n - 1, 0};
		HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
		                                  &coeffs[0][0]);
	}
	//middle
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
	// north
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		int stencil_indeces[3] = {0, 3, 4};
		vector<array<double, 3>> coeffs;
		array<double, 3>         this_coeffs;
		if (hasNbr(Side::north)) {
			this_coeffs[0] = -2.0 / (h_y * h_y);
			this_coeffs[1] = 1.0 / (h_y * h_y);
			this_coeffs[2] = 1.0 / (h_y * h_y);
		} else {
			this_coeffs[0] = -3.0 / (h_y * h_y);
			this_coeffs[1] = 1.0 / (h_y * h_y);
			this_coeffs[2] = 0;
		}
		for (int xi = 0; xi < n; xi++) {
			coeffs.push_back(this_coeffs);
		}
		int lower[2] = {0, n-1};
		int upper[2] = {n - 1, n-1};
		HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
		                                  &coeffs[0][0]);
	}
	// west
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		if (hasCoarseNbr(Side::west)) {
			if (isCoarseLeft(Side::west)) {
				int stencil_indeces[3] = {0, 1, 2};
				int off_stencil_indeces[3] = {5, 6, 7};
				vector<array<double, 3>> coeffs;
				vector<array<double, 3>> off_coeffs;
				array<double, 3>         this_coeffs_l;
				array<double, 3>         off_this_coeffs_l;
				this_coeffs_l[0] = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_l[1] = 0;
				this_coeffs_l[2] =4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_l[0] = 1.0 / 12.0 / (h_x * h_x);
				off_this_coeffs_l[1] = 1.0 / 2.0 / (h_x * h_x);
				off_this_coeffs_l[2] = -1.0 / 20.0 / (h_x * h_x);
				array<double, 3> this_coeffs_r;
				array<double, 3> off_this_coeffs_r;
				this_coeffs_r[0] = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_r[1] = 0;
				this_coeffs_r[2] = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_r[0] = -1.0 / 20.0 / (h_x * h_x);
				off_this_coeffs_r[1] = 1.0 / 2.0 / (h_x * h_x);
				off_this_coeffs_r[2] = 1.0 / 12.0 / (h_x * h_x);
				for (int yi = 0; yi < n/2; yi++) {
					coeffs.push_back(this_coeffs_l);
					coeffs.push_back(this_coeffs_r);
					off_coeffs.push_back(off_this_coeffs_l);
					off_coeffs.push_back(off_this_coeffs_r);
				}
				int lower[2] = {0, 0};
				int upper[2] = {0, n - 3};
				HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
				                                  &coeffs[0][0]);
				HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, off_stencil_indeces,
				                                  &off_coeffs[0][0]);
				array<double, 3> this_coeffs_uu;
				array<double, 3> off_this_coeffs_uu;
				this_coeffs_uu[0] = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_uu[1] = 0;
				this_coeffs_uu[2] = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_uu[0] = 1.0 / 12.0 / (h_x * h_x);
				off_this_coeffs_uu[1] = -3.0 / 10.0 / (h_x * h_x);
				off_this_coeffs_uu[2] = 3.0 / 4.0 / (h_x * h_x);
				int uu_idx[2]     = {0, n - 1};
				HYPRE_SStructMatrixAddToValues(A, ds.id, uu_idx, 0, 3, stencil_indeces,
				                               &this_coeffs_uu[0]);
				HYPRE_SStructMatrixAddToValues(A, ds.id, uu_idx, 0, 3, off_stencil_indeces,
				                               &off_this_coeffs_uu[0]);
				array<double, 3> this_coeffs_ul;
				array<double, 3> off_this_coeffs_ul;
				this_coeffs_ul[0] = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_ul[1] = 0;
				this_coeffs_ul[2] = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_ul[0] = -1.0 / 20.0 / (h_x * h_x);
				off_this_coeffs_ul[1] = 7.0 / 30.0 / (h_x * h_x);
				off_this_coeffs_ul[2] = 7.0 / 20.0 / (h_x * h_x);
				int ul_idx[2]     = {0, n - 2};
				HYPRE_SStructMatrixAddToValues(A, ds.id, ul_idx, 0, 3, stencil_indeces,
				                               &this_coeffs_ul[0]);
				HYPRE_SStructMatrixAddToValues(A, ds.id, ul_idx, 0, 3, off_stencil_indeces,
				                               &off_this_coeffs_ul[0]);

			} else {
				// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
				int stencil_indeces[3]     = {0, 1, 2};
				int off_stencil_indeces[3] = {5, 6, 7};
				vector<array<double, 3>> coeffs;
				vector<array<double, 3>> off_coeffs;
				array<double, 3>         this_coeffs_l;
				array<double, 3>         off_this_coeffs_l;
				this_coeffs_l[0]     = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_l[1]     = 0;
				this_coeffs_l[2]     = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_l[0] = 1.0 / 12.0 / (h_x * h_x);
				off_this_coeffs_l[1] = 1.0 / 2.0 / (h_x * h_x);
				off_this_coeffs_l[2] = -1.0 / 20.0 / (h_x * h_x);
				array<double, 3> this_coeffs_r;
				array<double, 3> off_this_coeffs_r;
				this_coeffs_r[0]     = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_r[1]     = 0;
				this_coeffs_r[2]     = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_r[0] = -1.0 / 20.0 / (h_x * h_x);
				off_this_coeffs_r[1] = 1.0 / 2.0 / (h_x * h_x);
				off_this_coeffs_r[2] = 1.0 / 12.0 / (h_x * h_x);
				for (int yi = 0; yi < n / 2; yi++) {
					coeffs.push_back(this_coeffs_l);
					coeffs.push_back(this_coeffs_r);
					off_coeffs.push_back(off_this_coeffs_l);
					off_coeffs.push_back(off_this_coeffs_r);
				}
				int lower[2] = {0, 2};
				int upper[2] = {0, n - 1};
				HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
				                                  &coeffs[0][0]);
				HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, off_stencil_indeces,
				                                  &off_coeffs[0][0]);
				array<double, 3> this_coeffs_lu;
				array<double, 3> off_this_coeffs_lu;
				this_coeffs_lu[0] = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_lu[1] = 0;
				this_coeffs_lu[2] = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_lu[0] = 7.0 / 20.0 / (h_x * h_x);
				off_this_coeffs_lu[1] = 7.0 / 30.0 / (h_x * h_x);
				off_this_coeffs_lu[2] = -1.0 / 20.0 / (h_x * h_x);
				int lu_idx[2]         = {0, 1};
				HYPRE_SStructMatrixAddToValues(A, ds.id, lu_idx, 0, 3, stencil_indeces,
				                               &this_coeffs_lu[0]);
				HYPRE_SStructMatrixAddToValues(A, ds.id, lu_idx, 0, 3, off_stencil_indeces,
				                               &off_this_coeffs_lu[0]);
				array<double, 3> this_coeffs_ll;
				array<double, 3> off_this_coeffs_ll;
				this_coeffs_ll[0]     = -4.0 / 3.0 / (h_x * h_x);
				this_coeffs_ll[1]     = 0;
				this_coeffs_ll[2]     = 4.0 / 5.0 / (h_x * h_x);
				off_this_coeffs_ll[0] = 3.0 / 4.0 / (h_x * h_x);
				off_this_coeffs_ll[1] = -3.0 / 10.0 / (h_x * h_x);
				off_this_coeffs_ll[2] = 1.0 / 12.0 / (h_x * h_x);
				int ll_idx[2]         = {0, 0};
				HYPRE_SStructMatrixAddToValues(A, ds.id, ll_idx, 0, 3, stencil_indeces,
				                               &this_coeffs_ll[0]);
				HYPRE_SStructMatrixAddToValues(A, ds.id, ll_idx, 0, 3, off_stencil_indeces,
				                               &off_this_coeffs_ll[0]);
			}
		} else {
			int stencil_indeces[3] = {0, 1, 2};
			vector<array<double, 3>> coeffs;
			array<double, 3>         this_coeffs;
			if (hasCoarseNbr(Side::west)) {
				this_coeffs[0] = -1.0 / (h_x * h_x);
				this_coeffs[1] = 0;
				this_coeffs[2] = 1.0 / (h_x * h_x);
			} else if (hasNbr(Side::west)) {
				this_coeffs[0] = -2.0 / (h_x * h_x);
				this_coeffs[1] = 1.0 / (h_x * h_x);
				this_coeffs[2] = 1.0 / (h_x * h_x);
			} else {
				this_coeffs[0] = -3.0 / (h_x * h_x);
				this_coeffs[1] = 0;
				this_coeffs[2] = 1.0 / (h_x * h_x);
			}
			for (int yi = 0; yi < n; yi++) {
				coeffs.push_back(this_coeffs);
			}
			int lower[2] = {0, 0};
			int upper[2] = {0, n - 1};
			HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
			                                  &coeffs[0][0]);
		}
	}
	// middle
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
	// east
	{
		// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
		if (hasFineNbr(Side::east)) {
			int stencil_indeces[5] = {0, 1, 2, 3, 4};
			int off_stencil_indeces[4] = {5, 6, 7, 8};
			vector<array<double, 5>> coeffs;
			vector<array<double, 4>> off_coeffs;
			array<double, 5> this_coeffs;
			array<double, 4> off_this_coeffs;
			this_coeffs[0] = -2.0 / (h_x * h_x);
			this_coeffs[1] = 1.0 / (h_x * h_x);
			this_coeffs[2] = 0;
			this_coeffs[3] = -1.0 / 30.0 / (h_x * h_x);
			this_coeffs[4] = -1.0 / 30.0 / (h_x * h_x);
			off_this_coeffs[0] = 1.0 / 3.0 / (h_x * h_x);
			off_this_coeffs[1] = 1.0 / 5.0 / (h_x * h_x);
			off_this_coeffs[2] = 1.0 / 3.0 / (h_x * h_x);
			off_this_coeffs[3] = 1.0 / 5.0 / (h_x * h_x);
			for (int yi = 0; yi < n; yi++) {
				coeffs.push_back(this_coeffs);
				off_coeffs.push_back(off_this_coeffs);
            }
			int lower[2] = {n - 1, 1};
			int upper[2] = {n - 1, n - 2};
			HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 5, stencil_indeces,
			                                  &coeffs[0][0]);
			HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 4, off_stencil_indeces,
			                                  &off_coeffs[0][0]);
			// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
			int off_stencil_indeces_ul[5] = {5, 6, 7, 8,9};
			array<double, 5> this_coeffs_u;
			array<double, 5> off_this_coeffs_u;
			this_coeffs_u[0] = -21.0/10.0 / (h_x * h_x);
			this_coeffs_u[1] = 1.0 / (h_x * h_x);
			this_coeffs_u[2] = 0;
			this_coeffs_u[3] = 1.0 / 15.0 / (h_x * h_x);
			this_coeffs_u[4] = 0;
			off_this_coeffs_u[0] = 1.0 / 3.0 / (h_x * h_x);
			off_this_coeffs_u[1] = 1.0 / 5.0 / (h_x * h_x);
			off_this_coeffs_u[2] = 1.0 / 3.0 / (h_x * h_x);
			off_this_coeffs_u[3] = 1.0 / 5.0 / (h_x * h_x);
			off_this_coeffs_u[4] = -1.0 / 30.0 / (h_x * h_x);
			int u_idx[2]         = {n - 1, n - 1};
			HYPRE_SStructMatrixAddToValues(A, ds.id, u_idx, 0, 5, stencil_indeces,
			                               &this_coeffs_u[0]);
			HYPRE_SStructMatrixAddToValues(A, ds.id, u_idx, 0, 5, off_stencil_indeces_ul,
			                               &off_this_coeffs_u[0]);
			array<double, 5> this_coeffs_l;
			array<double, 5> off_this_coeffs_l;
			this_coeffs_l[0] = -21.0/10.0 / (h_x * h_x);
			this_coeffs_l[1] = 1.0 / (h_x * h_x);
			this_coeffs_l[2] = 0;
			this_coeffs_l[3] = 0;
			this_coeffs_l[4] = 1.0 / 15.0 / (h_x * h_x);
			off_this_coeffs_l[0] = 1.0 / 3.0 / (h_x * h_x);
			off_this_coeffs_l[1] = 1.0 / 5.0 / (h_x * h_x);
			off_this_coeffs_l[2] = 1.0 / 3.0 / (h_x * h_x);
			off_this_coeffs_l[3] = 1.0 / 5.0 / (h_x * h_x);
			off_this_coeffs_l[4] = -1.0 / 30.0 / (h_x * h_x);
			int l_idx[2]         = {n - 1, 0};
			HYPRE_SStructMatrixAddToValues(A, ds.id, l_idx, 0, 5, stencil_indeces,
			                               &this_coeffs_l[0]);
			HYPRE_SStructMatrixAddToValues(A, ds.id, l_idx, 0, 5, off_stencil_indeces_ul,
			                               &off_this_coeffs_l[0]);
			/*int left_idx[2]    = {n - 1, n - 1};
			stencil_indeces[3] = 3;
			stencil_indeces[4] = 9;
			HYPRE_SStructMatrixAddToValues(A, ds.id, left_idx, 0, 9, stencil_indeces,
			                               &this_coeffs[0]);
			stencil_indeces[3] = 4;
			stencil_indeces[4] = 9;
			array<double, 9> this_coeffs_r;
			int stencil_indeces_r[8] = {0, 4, 9, 5, 6, 7, 8, 1};
			this_coeffs_r[0] = -11.0 / 10.0 / (h_x * h_x);
			this_coeffs_r[1] = -14.0/15.0 / (h_x * h_x);
			this_coeffs_r[2] = 0;
			this_coeffs_r[3] = -1.0 / 30.0 / (h_x * h_x);
			this_coeffs_r[4]         = 7.0 / 3.0 / (h_x * h_x);
			this_coeffs_r[5] = -4.0 / 5.0 / (h_x * h_x);
			this_coeffs_r[6] = 7.0 / 3.0 / (h_x * h_x);
			this_coeffs_r[7] = -4.0 / 5.0 / (h_x * h_x);
			this_coeffs_r[8] = 1.0 / (h_x * h_x);

			int right_idx[2]   = {n - 1, 0};
			HYPRE_SStructMatrixAddToValues(A, ds.id, right_idx, 0, 9, stencil_indeces_r,
			                               &this_coeffs_r[0]);
                                           */

		} else {
			int stencil_indeces[3] = {0, 1, 2};
			vector<array<double, 3>> coeffs;
			array<double, 3>         this_coeffs;
			if (hasFineNbr(Side::east)) {
				this_coeffs[0] = -1.0 / (h_x * h_x);
				this_coeffs[1] = 1.0 / (h_x * h_x);
				this_coeffs[2] = 0;
            }else
			if (hasNbr(Side::east)) {
				this_coeffs[0] = -2.0 / (h_x * h_x);
				this_coeffs[1] = 1.0 / (h_x * h_x);
				this_coeffs[2] = 1.0 / (h_x * h_x);
			} else {
				this_coeffs[0] = -3.0 / (h_x * h_x);
				this_coeffs[1] = 1.0 / (h_x * h_x);
				this_coeffs[2] = 0;
			}
			for (int yi = 0; yi < n; yi++) {
				coeffs.push_back(this_coeffs);
			}
			int lower[2] = {n - 1, 0};
			int upper[2] = {n - 1, n - 1};
			HYPRE_SStructMatrixAddToBoxValues(A, ds.id, lower, upper, 0, 3, stencil_indeces,
			                                  &coeffs[0][0]);
		}
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
}
