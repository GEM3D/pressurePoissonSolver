#include "Grid.h"

using namespace std;
Grid::Grid(double d_begin, double d_end, int m, double f(double x))
{
	this->m      = m;
	u            = LinkedVector<double>(m);
	uxx          = vector<double>(m);
	domain_begin = d_begin;
	domain_end   = d_end;
	h            = (d_end - d_begin) / m;

	for (int i = 0; i < m; i++) {
		double x  = domain_begin + (i + 0.5) / m * (domain_end - domain_begin);
		uxx.at(i) = f(x);
	}
}

double Grid::spaceDelta() { return h; }
int    Grid::size() { return m; }
bool   Grid::hasLeftNbr() { return u.left_nbr_ptr != nullptr; }
bool   Grid::hasRightNbr() { return u.right_nbr_ptr != nullptr; }
void Grid::setLeftNbr(Grid &nbr) { u.left_nbr_ptr = &nbr.u; }
void Grid::setRightNbr(Grid &nbr) { u.right_nbr_ptr = &nbr.u; }
