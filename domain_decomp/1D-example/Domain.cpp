#include <cmath>
#include "Domain.h"

using namespace std;
Domain::Domain(double d_begin, double d_end, int m, double f(double x))
{
	this->m      = m;
	u_curr       = valarray<double>(m);
	u_prev       = valarray<double>(m);
	u_xx         = valarray<double>(m);
	domain_begin = d_begin;
	domain_end   = d_end;
	h            = (d_end - d_begin) / m;

	for (int i = 0; i < m; i++) {
		double x = domain_begin + (i + 0.5) / m * (domain_end - domain_begin);
		u_xx[i]  = f(x);
	}
}

double Domain::spaceDelta() { return h; }
double Domain::domainBegin() { return domain_begin; }
double Domain::domainEnd() { return domain_end; }
int    Domain::size() { return m; }
bool   Domain::hasLeftNbr() { return left_gamma_ptr != nullptr; }
bool   Domain::hasRightNbr() { return right_gamma_ptr != nullptr; }
double Domain::rightGamma() { return *right_gamma_ptr; }
double Domain::leftGamma() { return *left_gamma_ptr; }
void   Domain::swapCurrPrev() { u_curr.swap(u_prev); }
