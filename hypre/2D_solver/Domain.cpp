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
	if (hasNbr(Side::north)) {
		int ilower[2]     = {0, n};
		int iupper[2]     = {n - 1, n};
		int nbr_ilower[2] = {0, 0};
		int nbr_iupper[2] = {n - 1, 0};
		int ret = HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::north),
		                                           nbr_ilower, nbr_iupper, index_map, index_dir);
		if (ret != 0) {
			throw std::runtime_error{"HYPRE_SStructGridSetNeighborPart returned: " + ret};
		}
	}
	if (hasNbr(Side::east)) {
		int ilower[2]     = {n, 0};
		int iupper[2]     = {n, n - 1};
		int nbr_ilower[2] = {0, 0};
		int nbr_iupper[2] = {0, n - 1};
		int ret = HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::east),
		                                           nbr_ilower, nbr_iupper, index_map, index_dir);
		if (ret != 0) {
			throw std::runtime_error{"HYPRE_SStructGridSetNeighborPart returned: " + ret};
		}
	}
	if (hasNbr(Side::south)) {
		int ilower[2]     = {0, -1};
		int iupper[2]     = {n - 1, -1};
		int nbr_ilower[2] = {0, n - 1};
		int nbr_iupper[2] = {n - 1, n - 1};
		int ret = HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::south),
		                                           nbr_ilower, nbr_iupper, index_map, index_dir);
		if (ret != 0) {
			throw std::runtime_error{"HYPRE_SStructGridSetNeighborPart returned: " + ret};
		}
	}
	if (hasNbr(Side::west)) {
		int ilower[2]     = {-1, 0};
		int iupper[2]     = {-1, n - 1};
		int nbr_ilower[2] = {n - 1, 0};
		int nbr_iupper[2] = {n - 1, n - 1};
		int ret = HYPRE_SStructGridSetNeighborPart(grid, ds.id, ilower, iupper, nbr(Side::west), nbr_ilower,
		                                 nbr_iupper, index_map, index_dir);
		if (ret != 0) {
			throw std::runtime_error{"HYPRE_SStructGridSetNeighborPart returned: " + ret};
        }
	}
}
void Domain::setMatrixCoeffs(HYPRE_SStructMatrix &A) {
	// int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	int stencil_indeces[5] = {0, 1, 2, 3, 4};
	array<double, 5> middle_coeffs = {-2.0 / (h_x * h_x) - 2.0 / (h_y * h_y), 1.0 / (h_x * h_x),
	                                  1.0 / (h_x * h_x), 1.0 / (h_y * h_y), 1.0 / (h_y * h_y)};
	array<double, 5> n_coeffs = {-2.0 / (h_x * h_x) - 3.0 / (h_y * h_y), 1.0 / (h_x * h_x),
	                             1.0 / (h_x * h_x), 1.0 / (h_y * h_y), 0 / (h_y * h_y)};
	array<double, 5> e_coeffs = {-3.0 / (h_x * h_x) - 2.0 / (h_y * h_y), 1.0 / (h_x * h_x),
	                             0 / (h_x * h_x), 1.0 / (h_y * h_y), 1.0 / (h_y * h_y)};
	array<double, 5> s_coeffs = {-2.0 / (h_x * h_x) - 3.0 / (h_y * h_y), 1.0 / (h_x * h_x),
	                             1.0 / (h_x * h_x), 0 / (h_y * h_y), 1.0 / (h_y * h_y)};
	array<double, 5> w_coeffs = {-3.0 / (h_x * h_x) - 2.0 / (h_y * h_y), 0 / (h_x * h_x),
	                             1.0 / (h_x * h_x), 1.0 / (h_y * h_y), 1.0 / (h_y * h_y)};
	array<double, 5> ne_coeffs = {-3.0 / (h_x * h_x) - 3.0 / (h_y * h_y), 1.0 / (h_x * h_x),
	                             1.0 / (h_x * h_x), 0 / (h_y * h_y), 0 / (h_y * h_y)};
	array<double, 5> se_coeffs = {-3.0 / (h_x * h_x) - 3.0 / (h_y * h_y), 1.0 / (h_x * h_x),
	                             0 / (h_x * h_x), 0 / (h_y * h_y), 1.0 / (h_y * h_y)};
	array<double, 5> sw_coeffs = {-3.0 / (h_x * h_x) - 3.0 / (h_y * h_y), 0 / (h_x * h_x),
	                             1.0 / (h_x * h_x), 0 / (h_y * h_y), 1.0 / (h_y * h_y)};
	array<double, 5> nw_coeffs = {-3.0 / (h_x * h_x) - 3.0 / (h_y * h_y), 0 / (h_x * h_x),
	                             1.0 / (h_x * h_x), 1.0 / (h_y * h_y), 0 / (h_y * h_y)};
	vector<array<double, 5>> coeffs;
	for (int yi = 0; yi < n; yi++) {
		for (int xi = 0; xi < n; xi++) {
			if (xi == n - 1 && yi == n - 1) {
				coeffs.push_back(ne_coeffs);
			} else if (xi == n - 1 && yi == 0) {
				coeffs.push_back(se_coeffs);
			} else if (xi == 0 && yi == 0) {
				coeffs.push_back(sw_coeffs);
			} else if (xi == 0 && yi == n-1) {
				coeffs.push_back(nw_coeffs);
            }else if(xi==0){
				coeffs.push_back(w_coeffs);
            }else if(xi==n-1){
				coeffs.push_back(e_coeffs);
            }else if(yi==0){
				coeffs.push_back(s_coeffs);
            }else if(yi==n-1){
				coeffs.push_back(n_coeffs);
			} else {
				coeffs.push_back(middle_coeffs);
			}
		}
	}
    int lower[2]={0,0};
    int upper[2]={n-1,n-1};
    cerr << coeffs.size() << endl;
	HYPRE_SStructMatrixSetBoxValues(A, ds.id, lower, upper, 0, 5, stencil_indeces, &coeffs[0][0]);
}
void Domain::fillRHS(HYPRE_SStructVector &b) {
	valarray<double> f_copy = f;
	if (!hasNbr(Side::north) && neumann) {
		f_copy[slice(n * (n - 1), n, 1)] -= 1 / h_y * boundary_north;
	} else {
		f_copy[slice(n * (n - 1), n, 1)] -= 2 / (h_y * h_y) * boundary_north;
	}
	if (!hasNbr(Side::east) && neumann) {
		f_copy[slice((n - 1), n, n)] -= 1 / h_x * boundary_east;
	} else {
		f_copy[slice((n - 1), n, n)] -= 2 / (h_x * h_x) * boundary_east;
	}
	if (!hasNbr(Side::south) && neumann) {
		f_copy[slice(0, n, 1)] += 1 / h_y * boundary_south;
	} else {
		f_copy[slice(0, n, 1)] -= 2 / (h_y * h_y) * boundary_south;
	}
	if (!hasNbr(Side::west) && neumann) {
		f_copy[slice(0, n, n)] += 1 / h_x * boundary_west;
	} else {
		f_copy[slice(0, n, n)] -= 2 / (h_x * h_x) * boundary_west;
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
void Domain::saveLHS(HYPRE_SStructVector &x)
{
	int lower[2] = {0, 0};
	int upper[2] = {n-1, n-1};
	HYPRE_SStructVectorGetBoxValues(x, ds.id, lower, upper, 0, &u[0]);
}
