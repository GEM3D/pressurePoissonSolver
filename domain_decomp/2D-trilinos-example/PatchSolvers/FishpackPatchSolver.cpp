#include "FishpackPatchSolver.h"
#include <valarray>
extern "C" {
void hstcrt_(double *a, double *b, int *m, int *mbdcnd, const double *bda, const double *bdb,
             double *c, double *d, int *n, int *nbdcnd, const double *bdc, const double *bdd,
             double *elmbda, double *f, int *idimf, double *pertrb, int *ierror, double *w);
}
using namespace std;
void FishpackPatchSolver::solve(Domain &d, const vector_type &f, vector_type &u,
                                const vector_type &gamma)
{
	double           h_x        = d.x_length / d.n;
	double           h_y        = d.y_length / d.n;
	auto             gamma_view = gamma.get1dView();
	auto             f_view     = f.get1dView();
	auto             u_view     = u.get1dViewNonConst();
	valarray<double> zeros(d.n);
	double           a      = d.x_start;
	double           b      = d.x_start + d.x_length;
	int              m      = d.n;
	int              mbcdnd = -1;
	if (d.isNeumann(Side::east) && d.isNeumann(Side::west)) {
		mbcdnd = 3;
	} else if (d.isNeumann(Side::west)) {
		mbcdnd = 4;
	} else if (d.isNeumann(Side::east)) {
		mbcdnd = 2;
	} else {
		mbcdnd = 1;
	}
	const double *bda = &zeros[0];
	const double *bdb = &zeros[0];

	double c      = d.y_start;
	double d2     = d.y_start + d.y_length;
	int    n      = d.n;
	int    nbcdnd = -1;
	if (d.isNeumann(Side::south) && d.isNeumann(Side::north)) {
		nbcdnd = 3;
	} else if (d.isNeumann(Side::south)) {
		nbcdnd = 4;
	} else if (d.isNeumann(Side::north)) {
		nbcdnd = 2;
	} else {
		nbcdnd = 1;
	}
	const double *bdc = &zeros[0];
	const double *bdd = &zeros[0];

	int     start  = d.id_local * d.n * d.n;
	double  elmbda = lambda;
	double *f_ptr  = &u_view[start];
	for (int i = 0; i < m * n; i++) {
		f_ptr[i] = f_view[start + i];
	}
	if (lambda != 0) {
		for (int i = 0; i < m * n; i++) {
			f_ptr[i] *=lambda;
		}
	}
	int              idimf  = n;
	if (d.hasNbr(Side::north)) {
		int idx = n * d.index(Side::north);
		for (int i = 0; i < n; i++) {
			f_ptr[n * (n - 1) + i] -= 2 / (h_y * h_y) * gamma_view[idx + i];
		}
	}
	if (d.hasNbr(Side::east)) {
		int idx = n * d.index(Side::east);
		for (int i = 0; i < n; i++) {
			f_ptr[i * n + (n - 1)] -= 2 / (h_x * h_x) * gamma_view[idx + i];
		}
	}
	if (d.hasNbr(Side::south)) {
		int idx = n * d.index(Side::south);
		for (int i = 0; i < n; i++) {
			f_ptr[i] -= 2 / (h_y * h_y) * gamma_view[idx + i];
		}
	}
	if (d.hasNbr(Side::west)) {
		int idx = n * d.index(Side::west);
		for (int i = 0; i < n; i++) {
			f_ptr[n * i] -= 2 / (h_x * h_x) * gamma_view[idx + i];
		}
	}
	double           pertrb = 0;
	int              ierror = 0;
	valarray<double> w(13 * m + 4 * n + m * log2(n));
	w[0] = 30.5;
	hstcrt_(&a, &b, &m, &mbcdnd, bda, bdb, &c, &d2, &n, &nbcdnd, bdc, bdd, &elmbda, f_ptr, &idimf,
	        &pertrb, &ierror, &w[0]);
	if (ierror != 0) {
		cerr << "IERROR: " << ierror << endl;
	}
}
