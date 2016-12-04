#include "Domain.h"

using namespace std;
Domain::Domain(double d_begin, double d_end, int m, double f(double x))
{
	this->m      = m;
	u            = vector<double>(m);
	uxx          = vector<double>(m);
	domain_begin = d_begin;
	domain_end   = d_end;
	h            = (d_end - d_begin) / m;

	for (int i = 0; i < m; i++) {
		double x  = domain_begin + (i + 0.5) / m * (domain_end - domain_begin);
		uxx.at(i) = f(x);
	}
}

double Domain::spaceDelta() { return h; }
int    Domain::size() { return m; }
bool   Domain::hasLeftNbr() { return left_nbr_ptr != nullptr; }
bool   Domain::hasRightNbr() { return right_nbr_ptr != nullptr; }
void Domain::setLeftNbr(Domain &nbr) { left_nbr_ptr = &nbr; }
void Domain::setRightNbr(Domain &nbr) { right_nbr_ptr = &nbr; }
vector<double> &Domain::getGrid(string str)
{
	if (str == "u_xx") {
		return uxx;
	}
	if (str == "u") {
		return u;
	}
}
vector<double> &Domain::getLeftGrid(string str) { return left_nbr_ptr->getGrid(str); }
vector<double> &Domain::getRightGrid(string str) { return right_nbr_ptr->getGrid(str); }
