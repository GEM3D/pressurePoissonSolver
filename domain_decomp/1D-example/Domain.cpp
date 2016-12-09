#include <cmath>
#include "Domain.h"

using namespace std;
Domain::Domain(double d_begin, double d_end, int m, double f(double x))
{
	this->m      = m;
	u_curr_ptr            = new vector<double>(m);
	u_prev_ptr            = new vector<double>(m);
	u_xx          = vector<double>(m);
	domain_begin = d_begin;
	domain_end   = d_end;
	h            = (d_end - d_begin) / m;

	for (int i = 0; i < m; i++) {
		double x  = domain_begin + (i + 0.5) / m * (domain_end - domain_begin);
		u_xx.at(i) = f(x);
	}
}

Domain::~Domain()
{
	delete u_curr_ptr;
	delete u_prev_ptr;
}

double Domain::spaceDelta() { return h; }
double Domain::domainBegin() { return domain_begin; }
double Domain::domainEnd() { return domain_end; }
int    Domain::size() { return m; }
bool   Domain::hasLeftNbr() { return left_nbr_ptr != nullptr; }
bool   Domain::hasRightNbr() { return right_nbr_ptr != nullptr; }
void Domain::setLeftNbr(Domain &nbr) { left_nbr_ptr = &nbr; }
void Domain::setRightNbr(Domain &nbr) { right_nbr_ptr = &nbr; }
vector<double> &                 Domain::uxx() { return u_xx; }
vector<double> &                 Domain::uCurr() { return *u_curr_ptr; }
vector<double> &                 Domain::uPrev() { return *u_prev_ptr; }
Domain &                         Domain::leftNbr() { return *left_nbr_ptr; }
Domain &                         Domain::rightNbr() { return *right_nbr_ptr; }
void                             Domain::swapCurrPrev()
{
	vector<double> *tmp = u_curr_ptr;
	u_curr_ptr          = u_prev_ptr;
	u_prev_ptr          = tmp;
}
double Domain::l2norm()
{
	double retval = 0;
	for (int i = 0; i < m; i++) {
		double diff = u_curr_ptr->at(i) - u_prev_ptr->at(i);
		retval += diff * diff;
	}
	return sqrt(retval);
}
